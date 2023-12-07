import argparse
import time
from Bio import SeqIO
from jit_open import Handle, Queue
from os.path import basename, splitext
from Bio.Seq import reverse_complement
import edlib

"""
Code for demultiplexing raw sequencing output for PTRE-seq library
"""

##################################################################################################################################################################
def match(fastq_file, barcode_handle, re_seq, file_format='fastq'):
	
	reader = SeqIO.index(fastq_file, file_format) # read fastq file

	read = []	
	read_type = []
	read_len = []
	reporters = []

	for record in reader:
		found = False

		seq = str(reader[record].seq)  

		read.append(reader[record].id)
		read_len.append(len(seq))
		
		best_barcode_dist = 1000
		best_re_dist = 1000
		best_reporter = ''
		best_read_type = ''

		for reporter, _barcode in barcodes.items():
			
			barcode = 'CGAG' + _barcode + 'GGTA' # add 4 nucleotide upstream and downstream of barcode  
			alignment = edlib.align(barcode, seq, mode = 'HW') # align barcode
			p = reporter.split('_')[0]	

			# if barcode perfectly matches, immediately assign read to that barcode/reporter
			if alignment['editDistance'] == 0:
				best_barcode_dist = alignment['editDistance']
				best_reporter = reporter
		
				reseq = re_seq[p]		
				alignment_re = edlib.align(reseq, seq, mode = 'HW')
				best_re_dist = alignment_re['editDistance']
				break		

			# if there is no perfect match, go through all barcode with 1 mismatch, align the 
			# regulatory element sequence to identify the most likely barcode	
			elif alignment['editDistance'] == 1:	
				# assign barcode dist to 1
				best_barcode_dist = 1
				reseq = re_seq[p]		
				alignment_re = edlib.align(reseq, seq, mode = 'HW')
				# loop to all barcode with 1 mismatch and find the RE with the best alignment
				if alignment_re['editDistance'] <= best_re_dist:
					best_re_dist = alignment_re['editDistance']
					best_reporter = reporter				

				found = True

		if best_barcode_dist == 0 and best_re_dist <= 5:
			best_read_type = 'match_zero'
		elif best_barcode_dist == 0 and best_re_dist > 5:
			best_read_type = 'mismatch_zero'
		elif best_barcode_dist == 1 and best_re_dist <= 5:
			best_read_type = 'match_one'
		elif best_barcode_dist == 1 and best_re_dist > 5:
			best_read_type = 'mismatch_one'
		else:
			best_read_type = 'cannot_find_barcode'			
			best_reporter = 'n/a'

		read_type.append(best_read_type)
		reporters.append(best_reporter)

	return read, read_type, read_len, reporters

################################################################################################

def parseArgs():
	parser = argparse.ArgumentParser(description='Demultiplex barcodes in fastq files')
	parser.add_argument('merged', type=str, help='Merged fastq file')
	parser.add_argument('barcode_file', type=str, help='Path to barcode file')	
	parser.add_argument('outfile', type=str, help='Output file')
	args = parser.parse_args()

	return args

#################################################################################################

if __name__== '__main__':
	
	start_time = time.time()
	args = parseArgs()

	# import sequence for regulatory elements
	all_re_seq = {}
	with open('all_re_seq.txt', 'r') as f:
		for line in f:
			spl = line.rstrip().split('\t')
			all_re_seq[spl[0]] = spl[1]	
	# import sequence for barcodes
	barcodes = {}
	with open(args.barcode_file, 'r') as f:
		for line in f:
			spl = line.rstrip().split('\t')
			barcodes[spl[0]] = spl[1]

	# align sequencing reads to barcodes and regulatory elements  
	read, read_type, read_len, reporters = match(args.merged, barcodes, all_re_seq)

	# write	output for fastq file
	with open(args.outfile, 'w') as out:
		for name, l, t, re in zip(read, read_len, read_type, reporters):
			out.write(f'{name}\t{l}\t{t}\t{re}\n')

	print(f'total time is {time.time()-start_time}')

	



























