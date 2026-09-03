# 13. Export of final sequences

**Objective:** extract the protein and CDS sequences of the final, non-redundant gene set of each
isolate.

**Output:** `${OUTDIR}/70_extracted_seq/`.

```bash
gffread -y ${OUTDIR}/70_extracted_seq/protein/${PREFIX}_protein.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ${OUTDIR}/60_merged/${PREFIX}_merged.gff

gffread -x ${OUTDIR}/70_extracted_seq/cds/${PREFIX}_cds.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ${OUTDIR}/60_merged/${PREFIX}_merged.gff
```

- `-y`: write protein sequences.
- `-x`: write CDS sequences.

The three files below constitute the final annotation of each isolate.

| File | Content |
|---|---|
| `${OUTDIR}/60_merged/${PREFIX}_merged.gff` | Non-redundant gene models with systematic identifiers |
| `${OUTDIR}/70_extracted_seq/protein/${PREFIX}_protein.fa` | Protein sequences of those gene models |
| `${OUTDIR}/70_extracted_seq/cds/${PREFIX}_cds.fa` | CDS sequences of those gene models |

The intermediate CDS files of [step 6](06_extract_cds.md) are compressed with `gzip` at the end of
the run.
