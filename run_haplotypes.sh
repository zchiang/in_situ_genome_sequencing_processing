# Author: Zachary Chiang, Buenrostro Lab, Harvard University
# Runs haplotyping pipeline for in situ genome sequencing data

# download VCF files

wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/C3H_HeJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz

# set file names

vcf_file_1='C3H_HeJ.mgp.v5.snps.dbSNP142.vcf.gz'
vcf_file_2='DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz'

bam_file='expt002_nova_all.mm10.high_pri.group.bam'
umi_file='expt002_nova_all.mm10.final.txt'

# get list of chromosomes

chrs=`samtools view -H $bam_file | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | awk '{if(length($0)<6)print}'`

# loop through chromosomes

for chr in $chrs; do

	echo $chr
	chr_num=${chr#"chr"}

	# finds haplotype-informative reads

	nohup python haplotype_reads_two_vcfs.py <(zcat $vcf_file_1 | grep -v -P "^#" | grep -P "^$chr_num\t") <(zcat $vcf_file_2 | grep -v -P "^#" | grep -P "^$chr_num\t")  <(samtools view -h $bam_file $chr) haplotyped_reads_nova/$chr\.bam &
	
	# sort resulting bam file
	
	nohup samtools sort -o haplotyped_reads_nova/$chr.st.bam haplotyped_reads_nova/$chr.bam &
#done

# merge all chromosomes

samtools merge -f exp002_nova_all.mm10.hap.bam haplotyped_reads_nova/*.st.bam

# add haplotype information to final bam file

python add_haplotypes.py $umi_file <(samtools view exp002_nova_all.mm10.hap.bam) > exp002_nova_all.mm10.final.haps.txt
