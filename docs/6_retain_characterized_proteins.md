
### 6. Retain characterized proteins (Helixer ↔ 70-15 RefSeq)

**Objective:** keep Helixer proteome models that show high similarity to characterized **70‑15** proteins.

```bash
# Extract Helixer proteome FASTA and make a DIAMOND DB
gffread -y results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff

diamond makedb --threads 1 \
  --in results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -d   results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa

# Search characterized 70‑15 proteins against that DB
diamond blastp --threads 1 --very-sensitive \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  -d results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -q 70-15_RefSeq_characterized_protein.fa \
  -o results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv

# Filter hits and recover matching Helixer IDs
./41_filter_diamond_results.py \
  results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv \
  results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  1> results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv \
  2> results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa

cut -f2 results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv | \
  sort | uniq \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt

gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff
```

**About `41_filter_diamond_results.py`:**

This script filters the DIAMOND BLASTP results to retain only Helixer gene models that are similar to characterized 70-15 proteins. It takes as input the DIAMOND tabular output and the Helixer proteome FASTA file. For each hit, it checks if the sequence identity is greater than 30% and the query (RefSeq) coverage is over 50%. If these criteria are met, the Helixer gene is considered a reliable match and is retained. The script outputs a tab-separated summary of matching pairs to stdout and the corresponding Heliker protein sequences to stderr. The resulting list of Helixer IDs is then used to extract the corresponding GFF features for downstream merging.
