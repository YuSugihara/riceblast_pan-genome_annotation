#! /usr/bin/env python3

import sys
import pyfastx


loci_file = open(sys.argv[1])
braker_proteome_file = pyfastx.Fasta(sys.argv[2])
helixer_secretome_file = pyfastx.Fasta(sys.argv[3])
miniprot_file = pyfastx.Fasta(sys.argv[4])

def get_seq(fasta, seq_id):
	return str(fasta[seq_id])

for line in loci_file:
	line = line.rstrip('\n')
	cols = line.split('\t')
	braker_ids = cols[3].split(',')
	helixer_ids = cols[4].split(',')
	miniprot_ids = cols[5].split(',')
	if miniprot_ids[0] != '-':
		for miniprot_id in miniprot_ids:
			already_exist = False
			match_braker_helixer = False
			miniprot_seq = get_seq(miniprot_file, miniprot_id)
			if braker_ids[0] != '-':
				for braker_id in braker_ids:
					braker_seq = get_seq(braker_proteome_file, braker_id)
					if braker_seq == miniprot_seq:
						already_exist = True
						break
					if helixer_ids[0] != '-':
						for helixer_id in helixer_ids:
							helixer_seq = get_seq(helixer_secretome_file, helixer_id)
							if helixer_seq == braker_seq:
								match_braker_helixer = True
								break
							if helixer_seq == miniprot_seq:
								already_exist = True
								break
			else:
				if helixer_ids[0] != '-':
					for helixer_id in helixer_ids:
						helixer_seq = get_seq(helixer_secretome_file, helixer_id)
						if helixer_seq == miniprot_seq:
							already_exist = True
							break
			if (not already_exist) and (not match_braker_helixer):
				if braker_ids[0] != '-' or helixer_ids[0] != '-':
					print(miniprot_id, line)
	