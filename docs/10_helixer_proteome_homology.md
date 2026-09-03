# 10. Homology screening of Helixer proteome models

**Objective:** restrict the non-overlapping Helixer proteome models of
[step 9](09_remove_overlapping_models.md) to those that match a characterized protein of
*M. oryzae* strain 70-15. Only this subset enters the merged annotation in
[step 12](12_merge_and_rename_ids.md).

**Output:** `${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff`.

```bash
# Extract the Helixer proteome sequences and build a DIAMOND database
gffread -y ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff

diamond makedb --threads 1 \
  --in ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -d   ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa

# Search the characterized 70-15 proteins against that database
diamond blastp --threads 1 --very-sensitive \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  -d ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -q data/70-15_RefSeq_characterized_protein.fa \
  -o ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv

# Filter the hits and collect the matching Helixer identifiers
scripts/41_filter_diamond_results.py \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  1> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv \
  2> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa

cut -f2 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv | \
  sort | uniq \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt

gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff
```

- `--very-sensitive`: DIAMOND search mode used to detect remote similarity.
- The characterized 70-15 proteins are used as the query and the Helixer models as the subject, so
  that `qcovhsp` reports the coverage of the characterized protein.

## `41_filter_diamond_results.py`

The script reads the tabular DIAMOND output and the Helixer proteome FASTA file. A hit is retained
when the sequence identity exceeds 30% and the query coverage, that is the coverage of the
characterized 70-15 protein, exceeds 50%. The Helixer models that satisfy these criteria
are taken to correspond to characterized proteins and are included in the merged annotation of
[step 12](12_merge_and_rename_ids.md).

The script writes a tab-separated summary of the retained pairs, comprising the RefSeq identifier,
the Helixer identifier, the query and subject coverage and the sequence identity, to standard
output, and the corresponding Helixer protein sequences in FASTA format to standard error. The
list of Helixer identifiers is then extracted from the summary with `cut -f2 | sort | uniq`.
