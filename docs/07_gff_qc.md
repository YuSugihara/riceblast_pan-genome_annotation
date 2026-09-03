# 7. Quality control of each GFF

**Objective:** tag the gene models with the structural problems they show, and record a summary
report per evidence track.

**Output:** `${OUTDIR}/20_gff_qc/${PREFIX}/`.

```bash
scripts/20_gff_qc.py <GFF> <CDS.fa> 150,180,195 10,25,50 \
  1> ${OUTDIR}/20_gff_qc/${PREFIX}/<source>_qc.gff \
  2> ${OUTDIR}/20_gff_qc/${PREFIX}/<source>_qc.txt
```

- Third argument: comma-separated CDS length thresholds in nucleotides (`150,180,195`).
- Fourth argument: comma-separated thresholds for the soft-masked fraction of the CDS, as a
  percentage (`10,25,50`).
- Standard output: the input GFF with the quality-control tags appended to each record.
- Standard error: a tabular summary of the tags assigned to each gene model.

The tags added by the script are the following.

| Tag | Meaning |
|---|---|
| `not_multiple_of_3` | The CDS length is not a multiple of three |
| `stop_codon_in_cds` | An in-frame stop codon occurs within the CDS |
| `no_start_codon` | The CDS does not begin with a start codon |
| `no_stop_codon` | The CDS does not end with a stop codon |
| `shorter_than_150nt` | The CDS is shorter than the given length threshold (one tag per threshold) |
| `masked_over_25` | The soft-masked fraction of the CDS exceeds the given threshold (one tag per threshold) |

The script is applied to each of the four evidence tracks produced in
[step 6](06_extract_cds.md); the tags are used as filtering criteria in
[step 8](08_filter_gene_models.md).
