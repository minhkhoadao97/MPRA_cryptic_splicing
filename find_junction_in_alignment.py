import pysam
from Bio import SeqIO
import pandas as pd
import sys, os
import re


def find_RT_error(cigar, forward_read):
	# function to identify RT error via complement nucleotide at position 289 in sequencing read
	
	cigar_tuple = re.findall(r'(\d+)([A-Z]|={1})', cigar)
	comp_found = False # check if there is a complementary sequence at position 289
	spacer_found = False # mark the position of the spacer sequence
	current_pos = 0
	
	for ele in reversed(cigar_tuple):
		if spacer_found:
			#print(cigar, current_pos, forward_read[-current_pos-3:-current_pos+6])
			if ele[1] == 'X' or ele[1] == 'D':
				if forward_read[-current_pos-3:-current_pos+6] == 'CGGATGCAT':
					comp_found = True
			break

		current_pos += int(ele[0])
		if ele[1] == '=' and forward_read[-current_pos:-current_pos+6] == 'ATGCAT':
			spacer_found = True

	if comp_found and spacer_found:
		return 'RT_error'
	else:
		return 'none'

def filter_cigar(cigar):
	# function to identify chimeric read via cigar string complexity

	cigar_dic = {'S': 0, 'D': 0, 'X': 0, 'M': 0, 'I': 0}
	cigar_p = re.findall(r'(\d+)([A-Z]{1})', cigar)
	for l, t in cigar_p:
		cigar_dic[t] += 1
	
	if cigar_dic['D'] >= 1 and cigar_dic['D'] + cigar_dic['X'] <= 5 and cigar_dic['I'] == 0:
		return 'pass'
	else:
		return 'fail'

def read_fasta(fn):
	seq = ''

	with open(fn, 'r') as f:
		for line in f:
			if not line.startswith('>'):
				seq += line.rstrip()

	return seq

bam_file = sys.argv[1]
bam = pysam.AlignmentFile(bam_file, 'r')

basename = bam_file.split('/')[-1]
sample = '_'.join(basename.split('.')[0].split('_')[0:1])
reporter = '_'.join(basename.split('.')[0].split('_')[-2:])
ref_seq = read_fasta(f'/storage/mustoe/home/mkdao/3UTR_L1_ref_seq/{reporter}.fa')
#ref_seq = read_fasta(f'{reporter}.fa')

df = pd.DataFrame(columns = ['read', 'cigar', 'RT_error', 'filtering', 'read_length', 'seq5', 'seq3', 'pos5', 'pos3', 'intron_length'])
df_idx = 0

for read in bam:
	if read.cigartuples is None:
		print(bam_file, ' is empty')
		exit()

	RT_error = find_RT_error(read.cigarstring, read.get_forward_sequence())
	#error_dict[read.query_name.split(' ')[0]] = RT_error
	#continue

	filtering = filter_cigar(read.cigarstring)

	current_pos = 0
	
	#read_name.append(read.query_name.split(' ')[0])
	_introns = []

	digits = re.findall('[0-9]+', read.cigarstring)
	chars = re.findall('[A-Z]+|=', read.cigarstring)
	tup = list(zip(chars, digits))

	del_found = False

	for idx, ele in enumerate(tup):

		if ele[0] == 'S':
			continue

		current_pos += int(ele[1])

		if ele[0] == 'D': # minimum intron length?
			start = current_pos - int(ele[1])
			end = current_pos
			
			if tup[idx-1][0] == '=' and  tup[idx+1][0] == '=' and int(ele[1]) > 4:
				del_found = True
				start_seq = ref_seq[start-6:start+2] # need to add 5 for the first 5 nucleotides				
				end_seq = ref_seq[end-2:end+6]
						
				# read name, cigar, RT error, pass/fail, read length, start seq, end seq, start pos, end pos, intron length	
				_row = [read.query_name.split(' ')[0], \
					read.cigarstring, \
					RT_error, \
					filtering, \
					len(read.get_forward_sequence()), \
					start_seq, \
					end_seq, \
					start, \
					end, \
					ele[1]]
				
				df.loc[df_idx] = _row
				df_idx += 1				
	if not del_found:
		_row = [read.query_name.split(' ')[0], \
			read.cigarstring, \
			RT_error, \
			filtering, \
			len(read.get_forward_sequence()), \
			'n/a', \
			'n/a', \
			'n/a', \
			'n/a', \
			'n/a']

		df.loc[df_idx] = _row
		df_idx += 1				


#df = pd.read_csv(f'JUNCTION_FROM_ALIGNMENT/{sample}_{reporter}_junction_sequences.txt', sep = '\t', \
#		names = ['read', 'cigar', 'RTerror', 'cigar_check', 'read_length', 'seq5', 'seq3', 'pos5', 'pos3', 'intron_length'])
#df['RTerror'] = df['read'].apply(lambda g: error_dict[g])
#print(df.sort_values(by='RTerror'))
df.to_csv(f'JUNCTION_FROM_ALIGNMENT/{sample}_{reporter}_junction_sequences.txt', sep = '\t', index = False, header = False) 
