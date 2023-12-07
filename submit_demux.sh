#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=2-00:00:00

for f in FASTQ/*fastq; do
	name=${f%.*}
	out=""$name"_demux.txt"
	sbatch demux_command.sh $f $out 
done
