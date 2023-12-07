#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --mem 1GB
#SBATCH --time=2-00:00:00

for f in ALIGNMENT/*; do
	python find_junction_in_alignment.py $f
done
