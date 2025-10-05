
## Running BRAKER

This section describes how to generate gene predictions for each genome assembly using [BRAKER](https://github.com/Gaius-Augustus/BRAKER).

### Table of contents:


### 1. Repeat masking

```bash
RepeatMasker -pa 4 -gff -nolow -xsmall -lib moryzae_repeats.fa.classified genome.fa
```

### 1. Repeat masking