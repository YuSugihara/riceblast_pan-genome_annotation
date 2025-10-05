### 7. Helixer‑unique secreted proteins overlapping BRAKER (50\_helixer\_uniq\_secreted\_proteins)

**Objective:** identify Helixer secretome models that **overlap** BRAKER loci (potential secreted proteins missing signal peptides in BRAKER proteins), then extract their sequences.

```bash
gffcompare -p ${PREFIX} -T \
  -o results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff
```

- `-p ${PREFIX}`: Sets the prefix for the output files.
- `-T`: Disables generation of `.tmap` and `.refmap` files.

```bash
./50_check_overlapping_effectors.py \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare.loci \
  ./Proteomes/${PREFIX}.proteins.fa \
  ./Secretomes/${PREFIX}_secretome.fa \
  ./Helixer/20_filtered/Proteomes/${PREFIX}_Helixer.proteome_filtered.fa \
  ./Helixer/20_filtered/Secretomes/${PREFIX}_Helixer.secretome_filtered.fa \
  1> results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.tsv \
  2> results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt
```

**About `50_check_overlapping_effectors.py`:**

This script identifies Helixer-predicted secreted proteins that overlap with BRAKER loci but are not present in the BRAKER secretome. 

```bash
# Select and translate those Helixer‑unique effectors
gffread --ids results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff

gffread -y results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.fa \
  -g ./Frozen_assemblies/${PREFIX}.fa \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff
```
