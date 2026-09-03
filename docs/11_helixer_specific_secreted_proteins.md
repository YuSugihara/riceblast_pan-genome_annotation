# 11. Helixer-specific secreted proteins at BRAKER loci

**Objective:** identify Helixer secretome models that overlap a BRAKER locus at which no BRAKER
model was predicted to be secreted. Such loci indicate genes for which the BRAKER model misses
the signal peptide, and the Helixer model is therefore preferred.

**Output:** `${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/`.

These models are not retained in [step 9](09_remove_overlapping_models.md), because they overlap a
BRAKER model and are removed by `bedtools subtract -A`.

```bash
gffcompare -p ${PREFIX} -T \
  -o ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff
```

- `-p ${PREFIX}`: prefix of the output files.
- `-T`: do not write the `.tmap` and `.refmap` files.
- The `.loci` file produced by gffcompare lists, for each locus, the BRAKER and Helixer
  transcripts that were assigned to it.

```bash
scripts/50_check_overlapping_effectors.py \
  ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare.loci \
  ./Proteomes/${PREFIX}.proteins.fa \
  ./Secretomes/${PREFIX}_secretome.fa \
  ./Helixer/20_filtered/Proteomes/${PREFIX}_Helixer.proteome_filtered.fa \
  ./Helixer/20_filtered/Secretomes/${PREFIX}_Helixer.secretome_filtered.fa \
  1> ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.tsv \
  2> ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt
```

## `50_check_overlapping_effectors.py`

The script takes the gffcompare `.loci` file together with four FASTA files: the BRAKER proteome,
the BRAKER secretome, the filtered Helixer proteome and the filtered Helixer secretome. For each
locus it retains the locus when both of the following conditions are met.

- None of the BRAKER transcripts of the locus is present in the BRAKER secretome, that is, no
  BRAKER model at this locus was predicted to be secreted.
- At least one Helixer transcript of the locus is present in the Helixer secretome.

The locus coordinates and the corresponding BRAKER and Helixer identifiers are written to standard
output as a tab-separated table, and the Helixer identifiers alone to standard error, which serves
as the input to `gffread --ids` below.

```bash
# Select and translate the Helixer-specific secreted proteins
gffread --ids ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff

gffread -y ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.fa \
  -g ./Frozen_assemblies/${PREFIX}.fa \
  ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff
```
