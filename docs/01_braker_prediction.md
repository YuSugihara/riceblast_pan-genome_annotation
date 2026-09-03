# 1. Gene prediction with BRAKER

**Objective:** soft-mask each genome assembly and predict gene models from RNA-seq and protein
evidence.

**Output:** `Frozen_assemblies/${PREFIX}.fa` (soft-masked assembly),
`Gene_Models/${PREFIX}.gff3` (gene models) and `Proteomes/${PREFIX}.proteins.fa` (their proteins).

## Contents

- [1.1. Repeat masking](#11-repeat-masking)
- [1.2. RNA-seq alignment](#12-rna-seq-alignment)
- [1.3. Gene prediction from RNA-seq evidence](#13-gene-prediction-from-rna-seq-evidence)
- [1.4. Gene prediction from protein evidence](#14-gene-prediction-from-protein-evidence)
- [1.5. Reconciliation of the two predictions](#15-reconciliation-of-the-two-predictions)
- [1.6. Identifier renaming and translation](#16-identifier-renaming-and-translation)

## 1.1. Repeat masking

Repeats were identified and soft-masked with RepeatMasker using a species-specific repeat library.

```bash
RepeatMasker -pa 4 -gff -nolow -xsmall -lib data/moryzae_repeats.fa.classified ${PREFIX}.fa
```

- `-pa 4`: number of parallel search jobs.
- `-gff`: write the annotation of repetitive elements in GFF format.
- `-nolow`: do not mask low-complexity regions or simple repeats.
- `-xsmall`: soft-mask repeats in lower case rather than replacing them with `N`.
- `-lib`: custom repeat library built from known *M. oryzae* repeat elements
  (`data/moryzae_repeats.fa.classified` in this repository).

The soft-masked assemblies are placed in `Frozen_assemblies/`. Soft masking is required by the
quality-control step, which reports the masked fraction of each CDS (see [step 7](07_gff_qc.md)).

## 1.2. RNA-seq alignment

Published RNA-seq reads were aligned to each assembly to provide intron hints for BRAKER. Four
libraries were used per isolate: two in planta time courses, one KE002 sample and one mycelium
sample.

```bash
hisat2-build -p 2 -f ${PREFIX}.fa ${PREFIX}

hisat2 -p 4 -x ${PREFIX} -1 ${R1} -2 ${R2} --max-intronlen 5000 \
  | samtools view -b -u -F 4 - \
  | samtools sort -o ${PREFIX}_${LIBRARY}_hisat2.sorted.bam -
samtools index ${PREFIX}_${LIBRARY}_hisat2.sorted.bam
```

- `--max-intronlen 5000`: maximum intron length, set for the compact fungal genome.
- `samtools view -F 4`: discard unmapped reads.

## 1.3. Gene prediction from RNA-seq evidence

```bash
braker.pl --fungus --species=magnaporthe_oryzae --gff3 --cores 4 --useexisting \
  --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
  --genome=${PREFIX}.fa \
  --bam=${PREFIX}_inplanta_A_hisat2.sorted.bam,${PREFIX}_inplanta_B_hisat2.sorted.bam,${PREFIX}_KE002_hisat2.sorted.bam,${PREFIX}_mycellium_hisat2.sorted.bam \
  --workingdir=braker1_${PREFIX}
```

- `--fungus`: use the GeneMark-ES branch point model for fungal genomes.
- `--species=magnaporthe_oryzae --useexisting`: reuse an existing AUGUSTUS parameter set rather
  than retraining for each isolate.
- `--bam`: the four alignment files of [1.2](#12-rna-seq-alignment).

## 1.4. Gene prediction from protein evidence

The same assemblies were annotated a second time from protein homology alone, using a query set
that combines the 70-15 reference proteome (UniProt proteome UP000009058) with a curated set of
known *M. oryzae* effectors and their homologues (`data/Mo_Effectors.fasta` in this repository).

```bash
braker.pl --fungus --species=magnaporthe_oryzae --gff3 --cores 4 --useexisting \
  --epmode --AUGUSTUS_ab_initio \
  --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
  --genome=${PREFIX}.fa \
  --prot_seq=70-15_MG8_eff_pep_newid.fa \
  --workingdir=braker2_${PREFIX}
```

- `--epmode`: run in protein-evidence mode with ProtHint.
- `--AUGUSTUS_ab_initio`: also produce the *ab initio* AUGUSTUS prediction alongside the
  hint-supported one.
- `--prot_seq`: the concatenation of the two sets above.

## 1.5. Reconciliation of the two predictions

The two gene sets were combined with TSEBRA using a configuration that gives precedence to the
protein-evidence run.

```bash
tsebra.py \
  -g braker1_${PREFIX}/augustus.hints.gtf,braker2_${PREFIX}/augustus.hints.gtf \
  -e braker1_${PREFIX}/hintsfile.gff,braker2_${PREFIX}/hintsfile.gff \
  -c data/prefer_braker2.cfg \
  -o ${PREFIX}_rnaseq+prot_hints_tsebra.gtf
```

`data/prefer_braker2.cfg` differs from the TSEBRA default in the hint-source weights and in
requiring a supported start and stop codon rather than supported introns:

| Parameter | Default | Used here | Meaning |
|---|---|---|---|
| `E` | 10 | 8 | weight of RNA-seq hints |
| `C` | 5 | 10 | weight of protein hints with a high spliced-alignment score |
| `P` | 0.1 | 0.1 | weight of protein hints with a lower spliced-alignment score |
| `intron_support` | 0.75 | 0 | required fraction of supported introns |
| `stasto_support` | 1 | 1 | required support for the start and stop codon |

## 1.6. Identifier renaming and translation

The reconciled GTF was renamed to carry the isolate prefix, converted to GFF3 and translated.

```bash
rename_gtf.py --gtf ${PREFIX}_rnaseq+prot_hints_tsebra.gtf --prefix ${PREFIX} \
  --out ${PREFIX}.renamed.gtf

gtf2gff.pl --gff3 --printExon --printIntron --includeStopInCDS \
  --out ${PREFIX}.gff3 < ${PREFIX}.renamed.gtf

gffread -g ${PREFIX}.fa -FSv -y ${PREFIX}.proteins.fa ${PREFIX}.gff3
```

- `rename_gtf.py` is distributed with TSEBRA; `gtf2gff.pl` with AUGUSTUS.
- `--includeStopInCDS`: keep the stop codon inside the CDS feature.
- `gffread -FSv -y`: write protein sequences, preserving the attributes of the input.

The resulting files are the `Gene_Models/${PREFIX}.gff3` and `Proteomes/${PREFIX}.proteins.fa`
used from [step 6](06_extract_cds.md) onwards. The secreted proteins among them are predicted in
[step 3](03_secretome_prediction.md).
