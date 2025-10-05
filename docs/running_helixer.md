
## Running Helixer

This section describes how to generate gene predictions for each genome assembly using [Helixer](https://github.com/weberlab-hhu/Helixer).

### Table of contents:
- [1. Prepare Input FASTA Files](#1-prepare-input-fasta-files)
- [2. Run Helixer with Singularity](#2-run-helixer-with-singularity)
- [3. Extract Protein and CDS Sequences](#3-extract-protein-and-cds-sequences)
- [4. Filter Helixer Annotations](#4-filter-helixer-annotations)
- [5. Output layout for Helixer annotation](#5-output-layout-for-helixer-annotation)

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