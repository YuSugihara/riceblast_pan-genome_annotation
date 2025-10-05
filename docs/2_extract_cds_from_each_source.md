### 2. Extract CDS from each source (10\_cds)

**Objective:** obtain CDS and protein FASTA per evidence track for downstream QC and cleanup. Especially, BRAKER generates gene models without start or stop codons, and non-triplet gene models.

```bash
gffread -x results/10_cds/${PREFIX}/${PREFIX}.braker.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Gene_Models/${PREFIX}.gff3

gffread -x results/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff

gffread -x results/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff

gffread -x results/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff
```

* `-x`: write CDS sequences
* `-y`: write protein sequences
* `-g`: reference genome for extraction.
