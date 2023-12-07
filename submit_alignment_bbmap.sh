#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=5-00:00:00

for f in DEMULTIPLEXED/*; do

	if [[ "$f" == *""* ]]; then
	
		
		reporter=$(echo ${f%.*} | cut -d '_' -f 3-5)
		ref_file=~/3UTR_L1_ref_seq/"$reporter".fa
		basename="$(basename $f .fastq)"
		out_sam=ALIGNMENT/"$basename".sam
	

		#while [ $(squeue -u mkdao | wc -l) -gt 11 ]; do
		#	sleep 1m
		#done

		sbatch alignment_command.sh $ref_file $f $out_sam $basename	
	fi
done
