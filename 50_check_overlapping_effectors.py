#! /usr/bin/env python3

import os
import sys
import pyfastx
import pandas as pd
from Bio.Seq import Seq


loci_file = pd.read_csv(sys.argv[1], sep='\t', header=None)
braker_proteome_fasta = pyfastx.Fasta(sys.argv[2])
braker_secretome_fasta = pyfastx.Fasta(sys.argv[3])
helixer_proteome_fasta =pyfastx.Fasta(sys.argv[4])
helixer_secretome_fasta =pyfastx.Fasta(sys.argv[5])

def get_seq(fasta, seq_id):
	if seq_id in fasta:
		return str(fasta[seq_id])

for i, row in loci_file.iterrows():
	braker_ids = row[3].split(',')
	helixer_id = row[4].split(',')
	if sum([1 if braker_id in braker_secretome_fasta else 0 for braker_id in braker_ids]) == 0:
		for braker_id in braker_ids:
			if braker_id != '-':
				for helixer_id in helixer_id:
					helixer_seq = get_seq(helixer_proteome_fasta, helixer_id)
					if helixer_id in helixer_secretome_fasta:
						print(row[1].replace('[+]', ':').replace('[-]', ':'), row[3], row[4], sep='\t', file=sys.stdout)
						print(row[4], file=sys.stderr)
						break