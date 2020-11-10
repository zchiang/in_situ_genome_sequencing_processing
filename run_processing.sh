# Author: Zachary Chiang, Buenrostro Lab, Harvard University
# Runs processing of ex situ sequencing data 

picard='/mnt/bin/picard/picard.jar'

genome=$1
prefix=$2
fastq_dir=$3
name=$4

genome_ver=`basename $genome`

# move to working directory and make subdirectories

cd /mnt/users/zack/projects/embryo_seq_expt002_nova

mkdir -p $fastq_dir
mkdir -p bams
mkdir -p tmp

# move UMI to R2

python bin/add_umis_from_R2.py $fastq_dir/$prefix\_R1*.fastq.gz $fastq_dir/$prefix\_R2*.fastq.gz $fastq_dir/$prefix\_R3*.fastq.gz $fastq_dir/$prefix\_R1.umi.fastq.gz $fastq_dir/$prefix\_R2.umi.fastq.gz

# trim adapters

python bin/pyadapter_trim.py -a $fastq_dir/$prefix\_R1.umi.fastq.gz -b $fastq_dir/$prefix\_R2.umi.fastq.gz
mv $prefix\_R*.umi.trim.fastq.gz $fastq_dir

# align
(bowtie2 -X2000 -p2 --rg-id $prefix -x $genome --very-sensitive -k5 -1 <(gunzip -c fastqs/$prefix\_R1.umi.trim.fastq.gz) -2 <(gunzip -c fastqs/$prefix\_R2.umi.trim.fastq.gz) \
      |	samtools view -bS - -o  bams/$name.$genome_ver.bam) 2>bams/$name.$genome_ver.align.log

# sort

java -jar -Djava.io.tmpdir=tmp $picard SortSam SO=coordinate I=bams/$name.$genome_ver.bam O=bams/$name.$genome_ver.sort.bam VALIDATION_STRINGENCY=SILENT

# move UMIs from header to tag

python bin/move_umis_to_UM_tag.py <(samtools view -h bams/$name.$genome_ver.sort.bam) | samtools view -bS - -o bams/$name.$genome_ver.rg.bam

# remove optical duplicates + estimate duplication rate

java -XX:ParallelGCThreads=16 -jar -Djava.io.tmpdir=/tmp/zack/ $picard MarkDuplicates INPUT=bams/$name.$genome_ver.rg.bam OUTPUT=bams/$name.$genome_ver.dedup.bam METRICS_FILE=bams/$name.$genome_ver.dups.log REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
samtools index bams/$name.$genome_ver.dedup.bam 

# split MAPQ >= 30 and MAPQ < 30
samtools view bams/$name.$genome_ver.dedup.bam -b -q 30 -U bams/$name.$genome_ver.low_all.bam > bams/$name.$genome_ver.high_all.bam
samtools view bams/$name.$genome_ver.high_all.bam -b -F 256 > bams/$name.$genome_ver.high_pri.bam
samtools index bams/$name.$genome_ver.high_pri.bam

# UMI-tools group

eval "$(/mnt/users/zack/anaconda3/bin/conda shell.bash hook)"
conda activate umi_tools_env # requires python 3
umi_tools group -I bams/$name.$genome_ver.high_pri.bam --extract-umi-method=tag --umi-tag=UM:Z --method=directional --edit-distance-threshold=2 --group-out=group/$name.$genome_ver.high_pri.group.tsv --log=group/$name.$genome_ver.high_pri.group.log --paired  --output-bam -S bams/$name.$genome_ver.high_pri.group.bam
conda deactivate

# get unique UMI groups and filter out Ns

cat group/$name.$genome_ver.high_pri.group.tsv | grep -v "contig" | cut -f2,3,7,8,9 | sort | uniq | sort -k5,5 -n | grep -P -v "N" > group/$name.$genome_ver.high.umis.txt

# filters

python bin/pos_filter.py < group/$name.$genome_ver.high.umis.txt | sort -k3,3 > group/$name.$genome_ver.filt1.txt
python bin/umi_filter.py group/$name.$genome_ver.filt1.txt group/$name.$genome_ver.filt2.txt group/$name.$genome_ver.multiumi_table.txt 
cat group/$name.$genome_ver.filt2.txt | sort -n -k5,5 > group/$name.$genome_ver.multiumi.txt

# add colorspace sequence to file

python bin/seq2colorspace_embryo_001.py < group/$name.$genome_ver.multiumi.txt > group/$name.$genome_ver.cs_reads.txt

# add empty haplotype information

python bin/empty_haplotype.py < group/$name.$genome_ver.cs_reads.txt > group/$name.$genome_ver.final.txt

exit

