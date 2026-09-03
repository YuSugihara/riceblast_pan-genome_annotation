# 5. miniprot searches

**Objective:** add gene models at loci left unannotated by the predictions of
[step 1](01_braker_prediction.md) and [step 2](02_helixer_prediction.md), by aligning the curated
secretome query set and the characterized 70-15 proteins to each genome assembly.

**Input:** `Frozen_assemblies/${PREFIX}.fa`, `all_secretomes.uniq.fa`,
`data/70-15_RefSeq_characterized_protein.fa` ([step 4](04_reference_protein_datasets.md)).

**Output:** `${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff`.

```bash
# Against the curated secretome panel
miniprot -t 2 -G 3k --gff -P SEC --outs=0.5 -p 0.3 \
  ./Frozen_assemblies/${PREFIX}.fa all_secretomes.uniq.fa \
  2> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "##PAF" > ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff

# Against the characterized 70-15 proteins
miniprot -t 2 -G 3k --gff -P RefSeq \
  ./Frozen_assemblies/${PREFIX}.fa data/70-15_RefSeq_characterized_protein.fa \
  2>> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "#" >> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff
```

- `-t 2`: number of threads per job.
- `-G 3k`: maximum intron length permitted.
- `--gff`: write the alignments in GFF3 format.
- `-P SEC` / `-P RefSeq`: prefix of the generated identifiers, which records the origin of the
  query protein and is carried through to the merged annotation.
- `--outs=0.5`: minimum ratio of the alignment score to the best score of the same query.
- `-p 0.3`: minimum ratio of the alignment score to the score of the best alignment in the same
  genomic region.
- The two searches are appended to a single GFF file; `grep -v` removes the PAF and header lines
  that miniprot interleaves with the GFF records.

The results of both searches are treated as a single miniprot evidence track in the subsequent
steps.
