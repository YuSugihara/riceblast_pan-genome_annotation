#! /usr/bin/env python

import sys


prefix = sys.argv[1]
gff = open(sys.argv[2])

gene_id = 0
for line in gff:
	if line.startswith('#'):
		if 'gffread' not in line:
			print(line, end='')
		continue
	line = line.rstrip('\n')
	cols = line.split('\t')
	feature = cols[2]
	if feature == 'locus':
		gene_id += 1
		transcript_id = 0
		cols[2] = 'gene'
		cols[8] = 'ID={0}_g{1}'.format(prefix, str(gene_id).zfill(5))
	elif feature == 'mRNA':
		transcript_id += 1
		exon_id = 1
		cds_id = 1
		cols[8] = 'ID={0}_g{1}.t{2};Parent={0}_g{1}'.format(prefix, str(gene_id).zfill(5), transcript_id)
	elif feature == 'exon':
		cols[8] = 'ID={0}_g{1}.t{2}.exon{3};Parent={0}_g{1}.t{2}'.format(prefix, str(gene_id).zfill(5), transcript_id, exon_id)
		exon_id += 1
	elif feature == 'CDS':
		cols[8] = 'ID={0}_g{1}.t{2}.CDS{3};Parent={0}_g{1}.t{2}'.format(prefix, str(gene_id).zfill(5), transcript_id, cds_id)
		cds_id += 1
	print('\t'.join(cols))
	previous_gene_id = gene_id
gff.close()