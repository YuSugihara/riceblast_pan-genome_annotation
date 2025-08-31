# Rice Blast Pan‑genome Annotation Methods Pipeline

[Pipeline overview (Figure 1)](figures/Figure_1_Annotation_pipeline.pdf)
*Figure 1. Overview of the end-to-end annotation pipeline (see Execution script and Pipeline workflow sections for commands).*

This repository provides a step-by-step annotation workflow and Materials & Methods (MM) for assembling a non-redundant gene set for each *Magnaporthe oryzae* genome. The pipeline integrates BRAKER, Helixer, and protein alignments (miniprot), producing a merged GFF3, CDS, and protein FASTA for each assembly. All commands and input/output layouts are included for reproducibility.

## Dependencies

The following software and Python packages are required to run the pipeline. Each program or package is listed on a separate line for clarity.

| Program / Package | Version  | Description / Link |
|-------------------|----------|--------------------|
| samtools          | 1.20     | [samtools](http://www.htslib.org/) |
| seqkit            | 2.8.1    | [seqkit](https://bioinf.shenwei.me/seqkit/) |
| gffread           | 0.12.7   | [gffread](https://github.com/gpertea/gffread) |
| miniprot          | 0.13     | [miniprot](https://github.com/lh3/miniprot) |
| bedtools          | 2.31.1   | [bedtools](https://bedtools.readthedocs.io/) |
| diamond           | 2.1.9    | [diamond](https://github.com/bbuchfink/diamond) |
| gffcompare        | 0.12.6   | [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) |
| Python 3          | 3.10.14  | [python.org](https://www.python.org/) |
| pandas            | 2.2.1    | [pandas](https://pandas.pydata.org/) |
| biopython         | 1.83     | [biopython](https://biopython.org/) |
| pyfastx           | 2.1.0    | [pyfastx](https://github.com/lmdu/pyfastx) |
| Helixer           |          | [Helixer](https://github.com/weberlab-hhu/Helixer) |
| Singularity       |          | [Singularity](https://sylabs.io/singularity/) |

---

## Table of Contents

**Sections under preparation by Joe Win:**
- [Repeat masking (placeholder)](#repeat-masking-placeholder)
- [BRAKER annotation (placeholder)](#braker-annotation-placeholder)
- [Signal peptide extraction (placeholder)](#signal-peptide-extraction-placeholder)

**Sections under preparation by Yu Sugihara:**
- [Running Helixer](#running-helixer)
  - [1. Prepare Input FASTA Files](#1-prepare-input-fasta-files)
  - [2. Run Helixer with Singularity](#2-run-helixer-with-singularity)
  - [3. Extract Protein and CDS Sequences](#3-extract-protein-and-cds-sequences)
  - [4. Filter Helixer Annotations](#4-filter-helixer-annotations)
  - [5. Output layout for Helixer annotation](#5-output-layout-for-helixer-annotation)
- [Dependencies](#dependencies)
- [Table of Contents](#table-of-contents)
- [Expected layout for inputs](#expected-layout-for-inputs)
- [Step-by-step pipeline](#step-by-step-pipeline)
  - [0. Build secretome dataset](#0-build-secretome-dataset)
    - [a. 70-15 _Magnaporthe oryzae_ secretome dataset](#a-70-15-magnaporthe-oryzae-secretome-dataset)
    - [b. Merge and deduplicate BRAKER and Helixer secretomes from all isolates](#b-merge-and-deduplicate-braker-and-helixer-secretomes-from-all-isolates)
    - [c. Filter 70‑15 RefSeq characterized proteins](#c-filter-70-15-refseq-characterized-proteins)
  - [1. Expected layout for outputs](#1-expected-layout-for-outputs)
  - [2. miniprot searches (00_miniprot)](#2-miniprot-searches-00_miniprot)
  - [3. Extract CDS from each source (10_cds)](#3-extract-cds-from-each-source-10_cds)
  - [4. QC each GFF (20_gff_qc)](#4-qc-each-gff-20_gff_qc)
  - [5. Filter gene models (30_filtered_gff)](#5-filter-gene-models-30_filtered_gff)
  - [6. Remove overlaps to define non-redundant gene sets (40_non_overlap_gff)](#6-remove-overlaps-to-define-non-redundant-gene-sets-40_non_overlap_gff)
  - [7. Retain conserved proteins (Helixer ↔ 70-15 RefSeq)](#7-retain-conserved-proteins-helixer--70-15-refseq)
  - [8. Helixer-unique effectors overlapping BRAKER (50_helixer_uniq_effectors)](#8-helixer-unique-effectors-overlapping-braker-50_helixer_uniq_effectors)
  - [9. Merge all evidence and rename IDs (60_merged)](#9-merge-all-evidence-and-rename-ids-60_merged)
  - [10. Export final sequences (70_extracted_seq)](#10-export-final-sequences-70_extracted_seq)
- [Parameter notes & tips](#parameter-notes--tips)
- [File naming conventions](#file-naming-conventions)
- [Appendix](#appendix)
  - [GFF and annotation file manipulation tips](#gff-and-annotation-file-manipulation-tips)

---

## Running Helixer

This section describes how to generate gene predictions for each genome assembly using [Helixer](https://github.com/weberlab-hhu/Helixer).

### 1. Prepare Input FASTA Files

Place all genome assemblies (with the `.fa` extension) in the `Frozen_assemblies/` directory.

### 2. Run Helixer with Singularity

Use the following script to process all assemblies in batch mode:

```bash
# Create required directories
mkdir -p ./Helixer/00_raw/GFF
mkdir -p ./Helixer/00_raw/Proteomes
mkdir -p ./Helixer/00_raw/Transcriptomes
mkdir -p ./Helixer/00_raw/Secretomes
mkdir -p ./Helixer/00_raw/Logs
mkdir -p ./Helixer/10_gff_qc
mkdir -p ./Helixer/20_filtered/GFF
mkdir -p ./Helixer/20_filtered/Proteomes
mkdir -p ./Helixer/20_filtered/Secretomes


FASTA_FILES=$(find ./Frozen_assemblies/*.fa)
for FASTA_FILE in $FASTA_FILES; do
  PREFIX=$(basename "$FASTA_FILE" .fa)
  singularity run --nv helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif \
    Helixer.py \
      --model-filepath fungi_v0.3_a_0300.h5 \
      --subsequence-length 21384 \
      --overlap-offset 10692 \
      --overlap-core-length 16038 \
      --species Helixer \
      --fasta-path "$FASTA_FILE" \
      --gff-output-path ./Helixer/00_raw/GFF/${PREFIX}_Helixer.gff \
      > ./Helixer/00_raw/Logs/${PREFIX}_Helixer.log
done
```

**Parameter notes:**
- `--model-filepath`: Use `fungi_v0.3_a_0300.h5` (best BUSCO performance for *Magnaporthe oryzae*; see [Figure S8](https://doi.org/10.1101/2023.02.06.527280))
- `--subsequence-length`: 21384 (recommended for fungi)
- `--overlap-offset`: 10692 (recommended for fungi)
- `--overlap-core-length`: 16038 (recommended for fungi)
- `--species`: Prefix for gene IDs
- `--fasta-path`: Input genome sequence
- `--gff-output-path`: Output GFF file

### 3. Extract Protein and CDS Sequences

Use `gffread` to extract protein and CDS FASTA files from the raw Helixer GFF:

```bash
function get_helixer_fasta() {
  PREFIX=$(basename "$1" .fa)
  gffread -y ./Helixer/00_raw/Proteomes/${PREFIX}_Helixer.protein.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ./Helixer/00_raw/GFF/${PREFIX}_Helixer.gff
  gffread -x ./Helixer/00_raw/Transcriptomes/${PREFIX}_Helixer.cds.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ./Helixer/00_raw/GFF/${PREFIX}_Helixer.gff
}
export -f get_helixer_fasta
find ./Frozen_assemblies/*.fa | xargs -P 4 -I"%" bash -c "get_helixer_fasta %"
```

### 4. Filter Helixer Annotations

Run the provided Python script to filter secretome and proteome sequences:

```bash
# Define a function to filter Helixer annotations
function filter_helixer_annotation() {
  PREFIX=$(basename "$1" .fa)
  ./00_filter_helixer_annotation.py ./Helixer/00_raw/GFF/${PREFIX}_Helixer.gff \
    ./Helixer/00_raw/Secretomes/${PREFIX}_Helixer_secretome.fa \
    2,3,4,6,11 \
    1> ./Helixer/10_gff_qc/${PREFIX}.helixer_pre_qc.gff \
    2> ./Helixer/10_gff_qc/${PREFIX}.helixer_pre_qc.txt

  # Filter secreted proteins and remove unwanted features
  grep 'secreted' ./Helixer/10_gff_qc/${PREFIX}.helixer_pre_qc.gff | \
    grep -v 'shorter_than_3bp' | grep -v 'five_prime_UTR' | grep -v 'three_prime_UTR' | grep -v 'exon' | \
    cut -f1-9 > ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff

  gffread -y ./Helixer/20_filtered/Secretomes/${PREFIX}_Helixer.secretome_filtered.fa \
    -g ./Frozen_assemblies/${PREFIX}.fa \
    ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff

  # Filter proteome (remove short exons and UTR/exon features)
  grep -v 'shorter_than_3bp' ./Helixer/10_gff_qc/${PREFIX}.helixer_pre_qc.gff | \
    grep -v 'five_prime_UTR' | grep -v 'three_prime_UTR' | grep -v 'exon' | \
    cut -f1-9 > ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff

  gffread -y ./Helixer/20_filtered/Proteomes/${PREFIX}_Helixer.proteome_filtered.fa \
    -g ./Frozen_assemblies/${PREFIX}.fa \
    ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff
}
export -f filter_helixer_annotation

# Run the filter function in parallel for all assemblies
find ./Frozen_assemblies/*.fa | xargs -P 4 -I"%" bash -c "filter_helixer_annotation %"
```

The script `00_filter_helixer_annotation.py` takes a GFF file and a secreted protein FASTA as input, and checks exon lengths (comma-separated list, e.g., `2,3,4,6,11`). It appends the following information to each GFF line:

```
<Transcript ID>\t<Transcript ID>.1\t<secreted>\t<exon shorter than threshold> ...
```

Filtered GFF and FASTA files for both secretome and proteome are then generated, with the following criteria:
- Genes with exons shorter than 3 bp are removed.
- Features such as `five_prime_UTR`, `three_prime_UTR`, and `exon` are excluded from the GFF.
- Only genes present in `${PREFIX}_Helixer_secretome.fa` are retained for the secretome.

### 5. Output layout for Helixer annotation

```
Helixer/00_raw/GFF/                 # Raw Helixer GFF (PREFIX_Helixer.gff)
Helixer/00_raw/Proteomes/           # Raw translated proteins (gffread)
Helixer/00_raw/Transcriptomes/      # Raw CDS (gffread)
Helixer/00_raw/Secretomes/          # Raw Helixer-predicted secretome FASTA (input to filter script)
Helixer/10_gff_qc/                  # QC‑annotated Helixer GFF + logs
Helixer/20_filtered/GFF/            # Filtered Helixer proteome & secretome GFF
Helixer/20_filtered/Proteomes/      # Filtered Helixer proteome FASTA
Helixer/20_filtered/Secretomes/     # Filtered Helixer secretome FASTA
```

---

## Expected layout for inputs

```
Frozen_assemblies/                  # Genome FASTA files; e.g., SAMPLE.fa
Gene_Models/                        # BRAKER or other curated GFF3 files; SAMPLE.gff3
Helixer/20_filtered/GFF/            # Filtered Helixer proteome & secretome GFF files
Helixer/20_filtered/Proteomes/      # Filtered Helixer proteome FASTA files
Helixer/20_filtered/Secretomes/     # Filtered Helixer secretome FASTA files
Proteomes/                          # BRAKER protein FASTA files for each sample
Secretomes/                         # BRAKER secretome FASTA files for each sample

effector_datasets/Yan_et_al_2023_dataset/
  ├── Yan_et_al_2023_Supplemental_Data_Set_S7.txt
  └── Magnaporthe_oryzae.MG8.pep.all.fa

effector_datasets/                  # Other curated effector/secretome panels
```


---

## Step-by-step pipeline

### 0. Build secretome dataset

#### a. 70-15 _Magnaporthe oryzae_ secretome dataset

```bash
cat Yan_et_al_2023_Supplemental_Data_Set_S7.txt | \
  samtools faidx -r - Magnaporthe_oryzae.MG8.pep.all.fa \
  > M_oryzae_MG8_XiaYan_secretome.fa
```

*This command reads protein IDs from the list and extracts the corresponding sequences from the 70-15 (MG8) proteome (Magnaporthe_oryzae.MG8.pep.all.fa).*

#### b. Merge and deduplicate BRAKER and Helixer secretomes from all isolates

```bash
cat Secretomes/*.fa Helixer/20_filtered/Secretomes/*.fa effector_datasets/*.fa > all_secretomes.fa
seqkit rmdup -s -i -P all_secretomes.fa > all_secretomes.uniq.fa
```

- `-s`: Remove duplicates based on sequence
- `-i`: Ignore case
- `-P`: Consider only the positive strand when comparing sequences

`effector_datasets/*.fa` contains secreted protein datasets generated by BRAKER.  
`Helixer/20_filtered/Secretomes/*.fa` contains secretome predictions from Helixer.

#### c. Filter 70‑15 RefSeq characterized proteins

```bash
./01_remove_uncharacterized_seq.py 70-15_RefSeq_protein.fa > 70-15_RefSeq_characterized_protein.fa
```

*Removes entries labelled “uncharacterized/hypothetical” to keep only characterized proteins for rescue.*


### 1. Expected layout for outputs

Per genome (where `PREFIX = basename(assembly .fa)`):

```
results/
  00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff                  # miniprot hits vs secretome & 70‑15 proteins
  10_cds/${PREFIX}/*.fa.gz                                      # CDS extracted per source (gzipped at end)
  20_gff_qc/${PREFIX}/*.txt, *.gff                              # QC reports and annotated GFFs
  30_filtered_gff/${PREFIX}/*.filtered.gff                      # Post‑QC filtered GFFs
  40_non_overlap_gff/${PREFIX}/*.{gff,txt,faa,tsv}              # Non‑overlap subsets; Helixer protein DB; DIAMOND results
  50_helixer_uniq_effectors/${PREFIX}/*                         # Helixer‑unique effectors & sequences
  60_merged/${PREFIX}_merged.gff                                # Final merged annotation (renamed IDs)
  70_extracted_seq/{protein,cds}/${PREFIX}_{protein,cds}.fa
```

---

### 2. miniprot searches (00\_miniprot)

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

### 3. Extract CDS from each source (10\_cds)

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

### 4. QC each GFF (20\_gff\_qc)

**Objective:** annotate and tabulate common issues (frameshifts, stop codons, short CDS, masked fraction) and keep a report.

```bash
./20_gff_qc.py <GFF> <CDS.fa> 150,180,195 10,25,50 \
  1> results/20_gff_qc/${PREFIX}/<source>_qc.gff \
  2> results/20_gff_qc/${PREFIX}/<source>_qc.txt
```

* Length thresholds (nt) to annotate: `150,180,195` (comma-separated list)
* Masked % thresholds to annotate: `10,25,50` (comma-separated list)
* The QC GFF adds tags like: `not_multiple_of_3`, `stop_codon_in_cds`, `no_start_codon`, `no_stop_codon`, `shorter_than_150nt`, `masked_over_25`.

### 5. Filter gene models (30\_filtered\_gff)

**Objective:** remove models that fail basic QC.

```bash
# BRAKER
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.braker_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff

# Helixer proteome
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_proteome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff

# Helixer secretome
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_secretome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff

# miniprot
grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.miniprot_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff
```

### 6. Remove overlaps to define non-redundant gene sets (40\_non\_overlap\_gff)

**Objective:** Construct a non-redundant set of gene models by prioritizing sources in the following order: **BRAKER** > **Helixer secretome** > **miniprot** > **Helixer proteome**.

- [`bedtools subtract`](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html) with `-A` removes entire features from the first input file (A) if they overlap with the second input file (B).
- `cut` is used to extract gene IDs from the GFF files.

1. **Helixer secretome vs BRAKER**

```bash
# Extract mRNA features from Helixer secretome
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff

# Remove Helixer secretome mRNAs overlapping BRAKER models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping Helixer secretome mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt
```

2. **miniprot vs BRAKER & Helixer secretome**

```bash
# Extract mRNA features from miniprot
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff

# Remove miniprot mRNAs overlapping BRAKER and Helixer secretome models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping miniprot mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt
```

3. **Helixer proteome vs BRAKER, Helixer secretome, and miniprot**

```bash
# Extract mRNA features from Helixer proteome
grep 'mRNA' results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff

# Remove Helixer proteome mRNAs overlapping BRAKER, Helixer secretome, and miniprot models
bedtools subtract -A \
  -a results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff \
  -b results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
bedtools subtract -A \
  -a stdin \
  -b results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff

# Extract IDs of non-overlapping Helixer proteome mRNAs
cut -f9 results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f1 | cut -d '=' -f2 \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt
```

4. **Extract non-redundant features**

- `gffread --ids` extracts features from a GFF file based on a list of IDs.

- `gffread -M -K --ids` extracts features while preserving the original GFF structure and attributes.

```bash
# Extract non-overlapping Helixer secretome GFF
gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff

# Extract non-overlapping miniprot GFF
gffread -M -K --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff

# Clean up miniprot annotation
./40_clean_up_miniprot_annotation.py \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff \
  results/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt

gffread -M -K --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff

# Extract non-overlapping Helixer proteome GFF
gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff
```

*Rationale:* Give precedence to BRAKER models, then add non-overlapping Helixer secretome models (likely effectors), then miniprot rescues, and finally additional Helixer proteome models.

### 7. Retain conserved proteins (Helixer ↔ 70-15 RefSeq)

**Objective:** keep Helixer proteome models that show high similarity to characterized **70‑15** proteins.

```bash
# Translate Helixer proteome subset and make a DIAMOND DB
gffread -y results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff

diamond makedb --threads 1 \
  --in results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -d   results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa

# Search characterized 70‑15 proteins against that DB
diamond blastp --threads 1 --very-sensitive \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  -d results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  -q 70-15_RefSeq_characterized_protein.fa \
  -o results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv

# Filter hits and recover matching Helixer IDs
./41_filter_diamond_results.py \
  results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv \
  results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
  1> results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv \
  2> results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa

cut -f2 results/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv | \
  sort | uniq \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt

gffread --ids results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff \
  > results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff
```

### 8. Helixer‑unique effectors overlapping BRAKER (50\_helixer\_uniq\_effectors)

**Objective:** identify Helixer secretome models that **overlap** BRAKER loci (potential effectors missing in BRAKER protein sets), then extract their sequences.

```bash
gffcompare -p ${PREFIX} -T \
  -o results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff

./50_check_overlapping_effectors.py \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare.loci \
  ./Proteomes/${PREFIX}.proteins.fa \
  ./Secretomes/${PREFIX}_secretome.fa \
  ./Helixer/20_filtered/Proteomes/${PREFIX}_Helixer.proteome_filtered.fa \
  ./Helixer/20_filtered/Secretomes/${PREFIX}_Helixer.secretome_filtered.fa \
  1> results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.tsv \
  2> results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt

# Select and translate those Helixer‑unique effectors
gffread --ids results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
  > results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff

gffread -y results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.fa \
  -g ./Frozen_assemblies/${PREFIX}.fa \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff
```

### 9. Merge all evidence and rename IDs (60_merged)

**Objective:** combine BRAKER, Helixer secretome non‑overlap, DIAMOND‑retained Helixer proteome, miniprot rescues, then give stable IDs.

```bash
# Pack Helixer evidence into one temp GFF
gffread -M -K \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff \
  > results/60_merged/temp/${PREFIX}_helixer_temp.gff

# Merge: BRAKER (base) + Helixer (temp) + miniprot cleaned
gffread --sort-alpha --cluster-only --force-exons \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  results/60_merged/temp/${PREFIX}_helixer_temp.gff \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff \
  > results/60_merged/temp/${PREFIX}_temp.gff

# Rename to stable per‑project IDs
./70_rename_ids.py ${PREFIX} results/60_merged/temp/${PREFIX}_temp.gff \
  > results/60_merged/${PREFIX}_merged.gff
```

### 10. Export final sequences (70\_extracted\_seq)

```bash
gffread -y results/70_extracted_seq/protein/${PREFIX}_protein.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/60_merged/${PREFIX}_merged.gff

gffread -x results/70_extracted_seq/cds/${PREFIX}_cds.fa \
        -g ./Frozen_assemblies/${PREFIX}.fa \
        results/60_merged/${PREFIX}_merged.gff

# Tidy: compress intermediate CDS sets from step 10
gzip results/10_cds/${PREFIX}/${PREFIX}.braker.fa \
     results/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa \
     results/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa \
     results/10_cds/${PREFIX}/${PREFIX}.miniprot.fa
```

---

## Appendix

### GFF and annotation file manipulation tips

This section provides practical tips and example commands for working with GFF3 and related annotation files, especially for extracting, filtering, and manipulating gene models.

#### Extracting features by attribute

- **Extract all mRNA features from a GFF3 file:**
  ```bash
  grep 'mRNA' input.gff > output.mrna.gff
  ```

#### Removing overlapping features

- **Remove features in file A that overlap any feature in file B using bedtools:**
  ```bash
  bedtools subtract -A -a fileA.gff -b fileB.gff > fileA.nonoverlap.gff
  ```
  - `-A`: Remove the entire feature in A if any overlap with B is found (not just the overlapping part).

#### Extracting IDs from GFF attributes

- **Extract the ID attribute from the 9th column of a GFF3 file:**
  ```bash
  cut -f9 file.gff | cut -d ';' -f1 | cut -d '=' -f2 > ids.txt
  ```

#### Extracting features by ID list

- **Extract features from a GFF file by a list of IDs using gffread:**
  ```bash
  gffread --ids ids.txt input.gff > output.filtered.gff
  ```
  - `--ids`: Provide a file with one ID per line.

- **Extract features and cluster transcripts into loci, removing redundant isoforms:**
  ```bash
  gffread -M -K --ids ids.txt input.gff > output.filtered.gff
  ```
  - `-M` / `--merge`: Cluster the input transcripts into loci, discarding "redundant" transcripts (those with the same exact introns and fully contained or equal boundaries).
  - `-K`: For use with `-M`, also discard as redundant the shorter, fully contained transcripts (intron chains matching a part of the container). This helps to keep only the longest or most representative isoform per locus.

#### Extracting sequences from GFF and genome

- **Extract CDS sequences:**
  ```bash
  gffread -x output.cds.fa -g genome.fa input.gff
  ```
  - `-x`: Output CDS FASTA.
  - `-g`: Reference genome.

- **Extract protein sequences:**
  ```bash
  gffread -y output.protein.fa -g genome.fa input.gff
  ```
  - `-y`: Output protein FASTA.
  - `-g`: Reference genome.

#### Deduplicating FASTA files

- **Remove duplicate sequences from a FASTA file using seqkit:**
  ```bash
  seqkit rmdup -s -i -P input.fa > output.uniq.fa
  ```
  - `-s`: Remove duplicates by sequence.
  - `-i`: Ignore case.
  - `-P`: Only consider the positive strand.

