### 9. Export final sequences (70\_extracted\_seq)

```bash
gffread -y results/70_extracted_seq/protein/${PREFIX}_protein.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/60_merged/${PREFIX}_merged.gff

gffread -x results/70_extracted_seq/cds/${PREFIX}_cds.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/60_merged/${PREFIX}_merged.gff
```
