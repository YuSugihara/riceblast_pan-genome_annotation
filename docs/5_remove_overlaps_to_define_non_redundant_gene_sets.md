
### 5. Remove overlaps to define non-redundant gene sets (40\_non\_overlap\_gff)

**Objective:** Construct a non-redundant set of gene models by prioritizing sources in the following order: **BRAKER** > **Helixer secretome** > **miniprot** > **Helixer proteome**.

### Table of contents:
- [a. Helixer secretome vs BRAKER](#a-helixer-secretome-vs-braker)
- [b. miniprot vs BRAKER & Helixer secretome](#b-miniprot-vs-braker--helixer-secretome)
- [c. Helixer proteome vs BRAKER, Helixer secretome, and miniprot](#c-helixer-proteome-vs-braker-helixer-secretome-and-miniprot)
- [d. Extract non-redundant features](#d-extract-non-redundant-features)

**Tools used:**
- [`bedtools subtract`](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html) with `-A` removes entire features from the first input file (A) if they overlap with the second input file (B).
- `cut` is used to extract gene IDs from the GFF files.

**a. Helixer secretome vs BRAKER**

```bash
# Extract mRNA features from Helixer secretome
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff

# Remove Helixer secretome mRNAs overlapping BRAKER models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping Helixer secretome mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt
```

**b. miniprot vs BRAKER & Helixer secretome**

```bash
# Extract mRNA features from miniprot
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff

# Remove miniprot mRNAs overlapping BRAKER and Helixer secretome models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping miniprot mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt
```

**c. Helixer proteome vs BRAKER, Helixer secretome, and miniprot**

```bash
# Extract mRNA features from Helixer proteome
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff

# Remove Helixer proteome mRNAs overlapping BRAKER, Helixer secretome, and miniprot models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping Helixer proteome mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt
```

**d. Extract non-redundant features**

- `gffread --ids` extracts features from a GFF file based on a list of IDs.

- `gffread -M -K --ids` extracts features while preserving the original GFF structure and attributes.

```bash
# Extract non-overlapping Helixer secretome GFF
gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff

# Extract non-overlapping miniprot GFF
gffread -M -K --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff

# Clean up miniprot annotation
./40_clean_up_miniprot_annotation.py \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff \
  results/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt

gffread -M -K --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff

# Extract non-overlapping Helixer proteome GFF
gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff
```

*Rationale:* Give precedence to BRAKER models, then add non-overlapping Helixer secretome models, then miniprot annotations, and finally additional Helixer proteome models.

**About `40_clean_up_miniprot_annotation.py`:**

This script is used to further filter the miniprot-derived gene models. It takes as input a GFF file (with miniprot annotations) and the corresponding CDS FASTA file. For each locus (typically a cluster of overlapping miniprot alignments), it checks whether the predicted CDS ends with a valid stop codon (`TAG`, `TAA`, or `TGA`). If at least one transcript in the locus has a valid stop codon, only those transcripts are retained; otherwise, all transcripts for that locus are kept. The script outputs a list of transcript IDs to retain, which is then used with `gffread --ids` to extract the cleaned-up features from the original GFF.