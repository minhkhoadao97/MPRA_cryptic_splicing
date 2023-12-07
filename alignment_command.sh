#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --mem=5000
#SBATCH --time=1-00:00:00

ref=$1
fastq=$2
out_sam=$3
basename=$4

~/softwares/bbmap/bbmap.sh ref=$ref in=$fastq out=$out_sam nodisk

out_bam=ALIGNMENT/"$basename".bam
sorted_bam=ALIGNMENT/"$basename".bam.sorted.bam
samtools view -S -b $out_sam > $out_bam
samtools sort $out_bam -o $sorted_bam
samtools index $sorted_bam
rm $out_sam $out_bam
