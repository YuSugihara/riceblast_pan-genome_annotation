#! /usr/bin/env python

import sys
import pyfastx


diamond_results = sys.argv[1]
helixer_fasta = pyfastx.Fasta(sys.argv[2])

seq_dict = {}
with open(diamond_results) as f:
	for line in f:
		line = line.rstrip()
		refseq_id = line.split('\t')[0]
		helixer_id = line.split('\t')[1]
		refseq_cov = float(line.split('\t')[12])/100
		helixer_cov = float(line.split('\t')[13])/100
		seq_idt = float(line.split('\t')[2])
		# If the sequence identity is over 30% and the RefSeq coverage is over 50%, the annotation is regarded as a reliable annotation.
		if seq_idt > 30 and refseq_cov > 0.5:
			if helixer_id not in seq_dict:
				seq_dict[helixer_id] = helixer_fasta[helixer_id].seq
			print(refseq_id, helixer_id, round(refseq_cov, 3), round(helixer_cov, 3), seq_idt, sep='\t', file=sys.stdout)
		
for k,v in seq_dict.items():
	print('>'+k, file=sys.stderr)
	print(v, file=sys.stderr)