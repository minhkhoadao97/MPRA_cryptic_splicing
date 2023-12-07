#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --mem=5000
#SBATCH --time=1-00:00:00

ref=$1 # reference sequence
fastq=$2 # demultiplexed fastq file
out_sam=$3 # name for sam file 
basename=$4 # name of reporter

# alignment command 
~/softwares/bbmap/bbmap.sh ref=$ref in=$fastq out=$out_sam nodisk

out_bam=ALIGNMENT/"$basename".bam
sorted_bam=ALIGNMENT/"$basename".bam.sorted.bam

# convert sam to bam files, sort and index
samtools view -S -b $out_sam > $out_bam
samtools sort $out_bam -o $sorted_bam
samtools index $sorted_bam
rm $out_sam $out_bam
