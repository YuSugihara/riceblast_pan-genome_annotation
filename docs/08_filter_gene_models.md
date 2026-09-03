# 8. Filtering of gene models

**Objective:** remove gene models that failed the quality control of
[step 7](07_gff_qc.md).

**Output:** `${OUTDIR}/30_filtered_gff/${PREFIX}/`.

The criteria differ between evidence tracks. In the table below, x marks a tag whose models were
removed from that track, and a dash a tag that was tolerated.

| Tag | BRAKER | Helixer proteome | Helixer secretome | miniprot |
|---|:--:|:--:|:--:|:--:|
| `not_multiple_of_3` | x | x | x | x |
| `stop_codon_in_cds` | x | x | x | x |
| `no_start_codon` | – | x | x | x |
| `no_stop_codon` | – | x | x | – |
| `shorter_than_150nt` | x | x | x | x |
| `masked_over_25` | x | x | x | x |

miniprot models that lack a stop codon are not removed here; they are resolved per locus in
[step 9](09_remove_overlapping_models.md), where transcripts terminating in a valid stop codon are
preferred whenever such a transcript exists at the locus.

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

`cut -f 1-9` removes the quality-control fields appended in [step 7](07_gff_qc.md) and restores a
valid nine-column GFF3 file.
