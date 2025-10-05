### 3. QC each GFF (20\_gff\_qc)

**Objective:** annotate and tabulate common issues (frameshifts, stop codons, short CDS, masked fraction) and keep a report.

```bash
./20_gff_qc.py <GFF> <CDS.fa> 150,180,195 10,25,50 \
  1> results/20_gff_qc/${PREFIX}/<source>_qc.gff \
  2> results/20_gff_qc/${PREFIX}/<source>_qc.txt
```

* Length thresholds (nt) to annotate: `150,180,195` (comma-separated list)
* Masked % thresholds to annotate: `10,25,50` (comma-separated list)
* The QC GFF adds tags like: `not_multiple_of_3`, `stop_codon_in_cds`, `no_start_codon`, `no_stop_codon`, `shorter_than_150nt`, `masked_over_25`.
