#! /usr/bin/env python3

import sys
import pyfastx
from Bio.Seq import Seq


miniprot_gff_file = open(sys.argv[1])
miniprot_fasta_file = pyfastx.Fasta(sys.argv[2])

def check_stop_codon(seq_id, miniprot_fasta_file=miniprot_fasta_file):
	cds_seq = str(miniprot_fasta_file[seq_id]).upper()
	if cds_seq[-3:] == 'TAG' or cds_seq[-3:] == 'TAA' or cds_seq[-3:] == 'TGA':
		return True
	else:
		return False


for line in miniprot_gff_file:
	if line.startswith('#'):
		continue
	line = line.rstrip('\n')
	cols = line.split('\t')
	if cols[1] == 'gffcl':
		seq_ids = cols[8].split('transcripts=')[-1].split(',')
		complete_seq_ids = [seq_id for seq_id in seq_ids if check_stop_codon(seq_id)]
		if len(complete_seq_ids) > 0:
			for seq_id in complete_seq_ids:
				print(seq_id)
		else:
			for seq_id in seq_ids:
				print(seq_id)
		
			

miniprot_gff_file.close()