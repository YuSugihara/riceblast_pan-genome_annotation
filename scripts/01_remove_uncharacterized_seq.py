#!/usr/bin/env python3

import sys
from Bio import SeqIO


fasta_file = sys.argv[1]

for record in SeqIO.parse(fasta_file, "fasta"):
	if "uncharacterized protein" not in record.description and \
	   "hypothetical protein" not in record.description:
		print(f">{record.description}")
		print(record.seq)