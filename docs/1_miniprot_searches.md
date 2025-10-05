
### 1. miniprot searches (00\_miniprot)

**Objective:** rescue models missed by ab initio sources by aligning curated secretomes and characterized proteins to the genome.

```bash
# Against curated secretome panel
miniprot -t 2 -G 3k --gff -P SEC --outs=0.5 -p 0.3 \
  $GENOME all_secretomes.uniq.fa \
  2> results/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "##PAF" > results/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff

# Against characterized 70‑15 proteins
miniprot -t 2 -G 3k --gff -P RefSeq \
  $GENOME 70-15_RefSeq_characterized_protein.fa \
  2>> results/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "#" >> results/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff
```

* `-t 2`: threads per job.
* `-G 3k`: max intron/gap length permitted.
* `--gff`: emit GFF3.
* `-P SEC|RefSeq`: labels the protein source; carried into attributes.
* `--outs=0.5 -p 0.3`: output and alignment score thresholds.