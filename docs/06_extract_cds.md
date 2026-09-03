# 6. CDS extraction from each source

**Objective:** extract the CDS sequences of the gene models of each evidence track, so that the
models can be evaluated by the quality-control step.

**Output:** `${OUTDIR}/10_cds/${PREFIX}/`.

The gene predictors can produce models that lack a start or a stop codon, or whose CDS length is
not a multiple of three. The extracted sequences are used in [step 7](07_gff_qc.md) to detect such
models.

```bash
gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.braker.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Gene_Models/${PREFIX}.gff3

gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff

gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff

gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff
```

- `-x`: write CDS sequences.
- `-g`: soft-masked genome assembly from which the sequences are extracted.

Four evidence tracks are produced: BRAKER, the filtered Helixer proteome, the filtered Helixer
secretome, and miniprot. Each track is treated separately in the quality control of
[step 7](07_gff_qc.md) and the filtering of [step 8](08_filter_gene_models.md).
