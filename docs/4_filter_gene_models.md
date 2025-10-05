
### 4. Filter gene models (30\_filtered\_gff)

**Objective:** remove models that fail basic QC.

```bash
# BRAKER
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.braker_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff

# Helixer proteome
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_proteome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff

# Helixer secretome
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_secretome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff

# miniprot
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.miniprot_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff
```