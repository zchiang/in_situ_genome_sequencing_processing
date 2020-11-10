# Author: Zachary Chiang, Buenrostro Lab, Harvard University
# Adds haplotype information to .bam file

import sys
import numpy as np
import math

group_file = open(sys.argv[1],'r')
sam_file = open(sys.argv[2],'r')

# create dict with all groups

groups = {}
for line in group_file:

	column = line.rstrip().split()
	column.append("-1")
	groups[int(column[4])] = column

# loop through bam file

sam_file = open(sys.argv[2],'r')

reads = {}
snp_counts = {}

for line in sam_file:

        # skip header lines

	if line[0] == "@":
                print line.rstrip()
                continue

        # paired-end read handling

        read = line.rstrip().split()
        read_name = read[0].split("_")[0]
        
	if read_name not in reads:
                reads[read_name] = line.rstrip()
        else:

                # set read info

                read_1 = reads[read_name].split()
                read_2 = read

                reads.pop(read[0],None)		

		umi_group = -1
		gts = []
	
                # look for haplotype tags

		for col in read_1:
			if "UG:i:" in col:
				umi_group = int(col.split(":")[2])
			elif "GT:Z:" in col:
				read_gts = ":".join(col.split(":")[2:]).split(";")
				gts = gts + read_gts
				
		for col in read_2:
			if "UG:i:" in col:
				umi_group = int(col.split(":")[2])
			elif "GT:Z:" in col:
				read_gts = ":".join(col.split(":")[2:]).split(";")
				gts = gts + read_gts

                frag_len = abs(int(read_1[8]))

                # make sure group is in final file 

                if umi_group in groups:
			
		        column = groups[umi_group]
	
                        # add info and remove duplicates

			if column[6] == "NA" and len(gts)>0:
				column[6] = set(gts)
                        elif len(gts) > 0:
				column[6] = set(list(column[6])+list(gts))

                        # increment observed reads count

                        for gt in gts:
                            if (umi_group,gt) in snp_counts:
                                snp_counts[(umi_group,gt)] += 1
                            else:
                                snp_counts[(umi_group,gt)] = 1
                       
                        # add frag length while we're here

                        if len(column) == 10:
                            column.append(str(frag_len))
                       
                        # add to dict for output

                        groups[umi_group] = column

# loop through output dict

for group,column in sorted(groups.items()):
	
        pass_snps = []
	
        # add info for haplotype-informative reads
        
        if column[6] != "NA":

                num_dups = int(column[3])
                umi_group = int(column[4])
                gts = column[6]

                # append all genotypes

                for gt in gts:
                    count = snp_counts[umi_group,gt]
                    if count > (float(num_dups)*0):
                        pass_snps.append(gt+":"+str(count))

                # add actual haplotype info

                if len(pass_snps) > 0:
                    

                    column[6] = ";".join(pass_snps)

                    # count maternal and paternal SNPs

                    for pass_snp in pass_snps:
                        hap = int(pass_snp.split(":")[2])
                        if hap == 0:
                            column[7] = str(int(column[7]) + 1)
                        elif hap == 1:
                            column[8] = str(int(column[8]) + 1)

                    # assign read to maternal, paternal, or conflicting

                    if int(column[7]) > 0 and int(column[8]) == 0:
                        column[9] = "0"
                    elif int(column[8]) > 0 and int(column[7]) == 0:
                        column[9] = "1"
                    else:
                        column[9] = "2"

	print "\t".join(column)



