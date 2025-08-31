#! /usr/bin/env python3

import os
import sys
import pandas as pd
from Bio import SeqIO


gff_name = sys.argv[1]
secretome_fasta_name = sys.argv[2]
exon_len_check_points = [int(n) for n in sys.argv[3].split(',')]

def read_gff(gff_name):
	gff = pd.read_csv(gff_name, sep='\t', header=None, comment='#', low_memory=False)
	gff = gff[gff[2]!='gene']
	gff['gene_id'] = gff[8].str.split(';').str[0].str.split('.').str[0].str.split('=').str[1]
	gff['transcript_id'] = gff[gff[2]!='mRNA'][8].str.split('Parent=').str[1].str.split(';').str[0]
	gff.loc[gff[2]=='mRNA', 'transcript_id'] = gff.loc[gff[2]=='mRNA'][8].str.split(';').str[0].str.replace('ID=', '')
	return gff

def read_fasta(protein_fasta_name):
	records_dict = SeqIO.to_dict(SeqIO.parse(protein_fasta_name, 'fasta'))
	return records_dict

def check_exon_len(gff):
	exon_len_dict = {}
	for i, row in gff.loc[gff[2]=='exon'].iterrows():
		exon_len = row[4] - row[3] + 1
		if row['transcript_id'] in exon_len_dict:
			exon_len_dict[row['transcript_id']].append(exon_len)
		else:
			exon_len_dict[row['transcript_id']] = [exon_len]
	return exon_len_dict

def check_gene_region(gff):
	gene_start_dict = {}
	gene_end_dict = {}
	for i, row in gff.loc[gff[2]=='CDS'].iterrows():
		start = row[3]
		end = row[4]
		if row['transcript_id'] in gene_start_dict:
			gene_start_dict[row['transcript_id']].append(start)
			gene_end_dict[row['transcript_id']].append(end)
		else:
			gene_start_dict[row['transcript_id']] = [start]
			gene_end_dict[row['transcript_id']] = [end]
	return gene_start_dict, gene_end_dict

def update_start_end(gff, gene_start_dict, gene_end_dict):
	for i, row in gff.loc[gff[2]=='mRNA'].iterrows():
		gff.loc[i, 3] = min(gene_start_dict[row['transcript_id']])
		gff.loc[i, 4] = max(gene_end_dict[row['transcript_id']])
	return gff

def check_qc(transcript_id, exon_len_dict, secretome_ids):
	qc_list = []
	if transcript_id in secretome_ids:
		qc_list.append('secreted')
	for exon_len_check_point in exon_len_check_points:
		if min(exon_len_dict[transcript_id]) < exon_len_check_point:
			qc_list.append('shorter_than_{}bp'.format(exon_len_check_point))
	return '|'.join(qc_list)


gff = read_gff(gff_name)
records_dict = read_fasta(secretome_fasta_name)
secretome_ids = [v.id for k, v in records_dict.items()]
exon_len_dict = check_exon_len(gff)
gene_start_dict, gene_end_dict = check_gene_region(gff)
gff = update_start_end(gff, gene_start_dict, gene_end_dict)

qc_lines = []
qc_lines_dict = {}
transcript_ids = list(gff['transcript_id'])
for transcript_id in transcript_ids:
	if transcript_id not in qc_lines_dict:
		qc_line = check_qc(transcript_id, exon_len_dict, secretome_ids)
		qc_lines_dict[transcript_id] = qc_line
		print(transcript_id, qc_line, sep='\t', file=sys.stderr)
	qc_lines.append(qc_lines_dict[transcript_id])

gff['qc'] = qc_lines
gff.to_csv(sys.stdout, sep='\t', header=None, index=False)