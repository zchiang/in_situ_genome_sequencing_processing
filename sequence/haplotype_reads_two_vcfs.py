# Author: Zachary Chiang, Buenrostro Lab, Harvard University
# Finds haplotype-informative reads 

import sys
import numpy as np
import os.path
import gzip
import pysam

snp_dict = {}
snp_chrs = [] 
snp_locs = [] 

# open files

vcf_file_1 = open(sys.argv[1],'r')
vcf_file_2 = open(sys.argv[2],'r')
bam_file = pysam.AlignmentFile(sys.argv[3],'rb')
out_file = pysam.AlignmentFile(sys.argv[4],'wb',template=bam_file)

# loop through VCFs, create dict (for hashed retrieval of ref/alt) and numpy array (for vectorized retrieval of SNPs in read)

for line in vcf_file_1:

    if line[0] == "#":
        continue

    column = line.rstrip().split()
    column[1] = str(int(column[1])-1) # VERY IMPORTANT, VCF IS 1-INDEXED, BAM IS 0-INDEXED

    snp_dict["chr"+column[0]+":"+column[1]] = column[3] + ":" + column[4] + ":" + column[9].split(":")[0] + ":1"
    snp_chrs.append("chr"+column[0])
    snp_locs.append(int(column[1]))

for line in vcf_file_2:

    if line[0] == "#":
        continue

    column = line.rstrip().split()
    column[1] = str(int(column[1])-1) # VERY IMPORTANT, VCF IS 1-INDEXED, BAM IS 0-INDEXED

    # CHECK IF SNP IS PRESENT IN BOTH

    if "chr"+column[0]+":"+column[1] in snp_dict:
        snp_dict["chr"+column[0]+":"+column[1]] = "both"
        continue

    snp_dict["chr"+column[0]+":"+column[1]] = column[3] + ":" + column[4] + ":" + column[9].split(":")[0] + ":2"
    snp_chrs.append("chr"+column[0])
    snp_locs.append(int(column[1]))

# convert to np arrays

snp_chrs = np.asarray(snp_chrs).astype(str)
snp_locs = np.asarray(snp_locs)

# loop through bam file

reads = {}

for read in bam_file.fetch():

    read_name = read.query_name

    # paired-end read handling

    if read_name not in reads:
	reads[read_name] = read
        continue
    else:
    	r1 = reads.pop(read_name)
	r2 = read
  
    # find overlapping haplotype-informative SNPs

    r1_snps = np.flatnonzero(np.logical_and(snp_chrs==r1.reference_name,np.logical_and(snp_locs >= r1.reference_start, snp_locs < r1.reference_start+r1.query_length-1)))
    r2_snps = np.flatnonzero(np.logical_and(snp_chrs==r2.reference_name,np.logical_and(snp_locs >= r2.reference_start, snp_locs < r2.reference_start+r1.query_length-1)))

    snps = []

    # loop through read 1 SNPs

    for i in range(len(r1_snps)):

	chrom = snp_chrs[r1_snps[i]]
	pos = snp_locs[r1_snps[i]]
        
        # SNP is in both VCF, likely not informative

        if snp_dict[chrom+":"+str(pos)] == "both":
            continue

        # get actual base

        gt = snp_dict[chrom+":"+str(pos)].split(":")
	base = r1.query_sequence[pos-r1.reference_start]
        
        # base is not either ref or alt

        if base != gt[0] and base != gt[1]:
            continue
        
        # base is ref, not helpful for non-reference F1 crosses

        if base == gt[0]:
            continue

        # base is alt, figure out which parent it matches

        elif base == gt[1] and (gt[3]) == "1":
            hap = 0
        elif base == gt[1] and gt[3] == "2":
            hap = 1

        # append to list of SNPs

        snps.append(str(pos)+":"+base+":"+str(hap))

    # same as above but for read 2

    for i in range(len(r2_snps)):

	chrom = snp_chrs[r2_snps[i]]
	pos = snp_locs[r2_snps[i]]
        
        if snp_dict[chrom+":"+str(pos)] == "both":
            break
        
        gt = snp_dict[chrom+":"+str(pos)].split(":")
	base = r2.query_sequence[pos-r2.reference_start]
        
        hap = (base == gt[1])
        if base != gt[0] and base != gt[1]:
            continue
        
        if base == gt[0]:
            continue
        elif base == gt[1] and gt[3] == "1":
            hap = 0
        elif base == gt[1] and gt[3] == "2":
            hap = 1

        snps.append(str(pos)+":"+base+":"+str(hap))
   
    # deduplicate SNPs (if insert is small a SNP will appear on both reads)

    snps =  np.unique(snps)

    # add SNPs to GT tag

    if len(snps)>0:
        r1.set_tag('GT',';'.join(snps),'Z')
        r2.set_tag('GT',';'.join(snps),'Z')

    # write reads to output

    out_file.write(r1)
    out_file.write(r2)
    
