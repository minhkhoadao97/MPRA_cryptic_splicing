import pandas as pd
from Bio import SeqIO
import sys, os

report_file = sys.argv[1]
name = report_file.split('.')[0]
fastq_file = sys.argv[2]

df = pd.read_csv(report_file, sep = '\t', usecols = [0, 3], names = ['read_name', 'best_reporter'])
reader = SeqIO.index(fastq_file, 'fastq')

for reporter in df['best_reporter'].unique():
	temp = df[df['best_reporter']==reporter]

	with open(f'DEMULTIPLEXED/{name}_{reporter}.fastq', 'a') as handle:
		for read in temp['read_name']:
			SeqIO.write(reader[read], handle, 'fastq')

