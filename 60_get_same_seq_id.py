#! /usr/bin/env python

import sys
import pyfastx


query_fasta_file = pyfastx.Fasta(sys.argv[1])
subject_fasta_file = pyfastx.Fasta(sys.argv[2])

query_seqs = set([record.seq for record in query_fasta_file])

for record in subject_fasta_file:
	if record.name.startswith('SEC'):
		if record.seq in query_seqs:
			print(record.name, sep='\n', file=sys.stdout)
