#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --mem=1
#SBATCH --time=2-00:00:00

# 1 = merged fastq file 
# 2 = output txt file 
python demux_barcode_re.py $1 ~/barcode_indexed.csv $2 
