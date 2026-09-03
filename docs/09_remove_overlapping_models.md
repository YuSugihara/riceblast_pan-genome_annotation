# 9. Removal of overlapping gene models

**Objective:** construct a non-redundant set of gene models by resolving overlaps between the
four evidence tracks according to a fixed source priority.

**Output:** `${OUTDIR}/40_non_overlap_gff/${PREFIX}/`.

The source priority is **BRAKER > Helixer secretome > miniprot > Helixer proteome**. BRAKER models
are retained in full; models of each subsequent track are retained only where they do not overlap
any model of a higher-priority track. Predicted secreted proteins are given precedence over the
remaining Helixer models because they are the focus of the study.

Passing this step is not sufficient for the Helixer proteome track. Its non-overlapping models are
filtered further in [step 10](10_helixer_proteome_homology.md), and only those matching a
characterized 70-15 protein enter the merged annotation.

Two tools are used throughout this step.

- [`bedtools subtract -A`](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html)
  removes a feature of file A in its entirety whenever it overlaps any feature of file B.
- `cut` extracts the identifiers from the attribute column of the resulting GFF files.

## Contents

- [9.1. Helixer secretome against BRAKER](#91-helixer-secretome-against-braker)
- [9.2. miniprot against BRAKER and the Helixer secretome](#92-miniprot-against-braker-and-the-helixer-secretome)
- [9.3. Helixer proteome against the three higher-priority tracks](#93-helixer-proteome-against-the-three-higher-priority-tracks)
- [9.4. Extraction of the non-redundant features](#94-extraction-of-the-non-redundant-features)

## 9.1. Helixer secretome against BRAKER

```bash
# Extract mRNA features from the Helixer secretome
grep 'mRNA' ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff

# Remove Helixer secretome mRNAs overlapping BRAKER models
bedtools subtract -A \
  -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff \
  -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff

# Extract the identifiers of the non-overlapping Helixer secretome mRNAs
cut -f9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt
```

## 9.2. miniprot against BRAKER and the Helixer secretome

```bash
# Extract mRNA features from the miniprot alignments
grep 'mRNA' ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff

# Remove miniprot mRNAs overlapping BRAKER or Helixer secretome models
bedtools subtract -A \
  -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff \
  -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff

# Extract the identifiers of the non-overlapping miniprot mRNAs
cut -f9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt
```

## 9.3. Helixer proteome against the three higher-priority tracks

```bash
# Extract mRNA features from the Helixer proteome
grep 'mRNA' ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff

# Remove Helixer proteome mRNAs overlapping BRAKER, Helixer secretome or miniprot models
bedtools subtract -A \
  -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff \
  -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
bedtools subtract -A \
  -a stdin \
  -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff

# Extract the identifiers of the non-overlapping Helixer proteome mRNAs
cut -f9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt
```

## 9.4. Extraction of the non-redundant features

The identifier lists produced above are used to extract the corresponding features from the
filtered GFF files of [step 8](08_filter_gene_models.md).

- `gffread --ids` extracts the features listed in a file of identifiers.
- `gffread -M -K --ids` additionally clusters the transcripts into loci and discards redundant
  transcripts, which is necessary for the miniprot track because several query proteins can align
  to the same locus.

```bash
# Non-overlapping Helixer secretome models
gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff

# Non-overlapping miniprot models
gffread -M -K --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff

# Select miniprot transcripts terminating in a valid stop codon
scripts/40_clean_up_miniprot_annotation.py \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff \
  ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt

gffread -M -K --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff

# Non-overlapping Helixer proteome models
gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff
```

### `40_clean_up_miniprot_annotation.py`

The script takes the non-overlapping miniprot GFF file and the corresponding CDS FASTA file of
[step 6](06_extract_cds.md). For each locus, which typically consists of a cluster of overlapping
miniprot alignments, it tests whether the predicted CDS terminates in a valid stop codon (`TAG`,
`TAA` or `TGA`). If at least one transcript of the locus satisfies this condition, only those
transcripts are retained; otherwise all transcripts of the locus are retained. The retained
transcript identifiers are written to standard output and are used with `gffread --ids` to
extract the cleaned-up features from the filtered GFF file.
