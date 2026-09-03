# 2. Gene prediction with Helixer

**Objective:** generate gene models for each genome assembly with
[Helixer](https://github.com/weberlab-hhu/Helixer) and derive filtered proteome and secretome
sets from them.

**Output:** `Helixer/` (see [2.5](#25-output-layout)).

## Contents

- [2.1. Input files](#21-input-files)
- [2.2. Running Helixer](#22-running-helixer)
- [2.3. Extraction of protein and CDS sequences](#23-extraction-of-protein-and-cds-sequences)
- [2.4. Filtering of the Helixer annotation](#24-filtering-of-the-helixer-annotation)
- [2.5. Output layout](#25-output-layout)

## 2.1. Input files

All soft-masked genome assemblies are placed in `Frozen_assemblies/` with the `.fa` extension; they
are produced by the repeat masking of [step 1](01_braker_prediction.md#11-repeat-masking).
The filtering script of [2.4](#24-filtering-of-the-helixer-annotation) additionally requires the
secreted proteins of each isolate in `Helixer/00_raw/Secretomes/${PREFIX}_Helixer_secretome.fa`;
these are predicted from the proteomes of [2.3](#23-extraction-of-protein-and-cds-sequences) as
described in [step 3](03_secretome_prediction.md).

## 2.2. Running Helixer

Helixer was run for all assemblies through its Singularity image.

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

- `--model-filepath`: `fungi_v0.3_a_0300.h5`, the fungal model that achieved the highest BUSCO
  completeness for *M. oryzae* (see Figure S8 of the [Helixer benchmarking study](https://doi.org/10.1101/2023.02.06.527280)).
- `--subsequence-length`, `--overlap-offset`, `--overlap-core-length`: values recommended for
  fungal genomes.
- `--species`: prefix used for the predicted gene identifiers.
- `--fasta-path`: input genome assembly.
- `--gff-output-path`: output GFF file.

## 2.3. Extraction of protein and CDS sequences

Protein and CDS sequences were extracted from the raw Helixer GFF with `gffread`.

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

- `-y`: write protein sequences.
- `-x`: write CDS sequences.
- `-g`: reference genome from which the sequences are extracted.

## 2.4. Filtering of the Helixer annotation

`00_filter_helixer_annotation.py` takes a Helixer GFF file and the corresponding secreted protein
FASTA file, and appends the following fields to each line of the GFF:

```
<Transcript ID>\t<Transcript ID>.1\t<secreted>\t<exon shorter than threshold> ...
```

The third argument is a comma-separated list of exon-length thresholds in base pairs. The
annotated GFF is then filtered into a proteome set and a secretome set.

```bash
# Define a function to filter Helixer annotations
function filter_helixer_annotation() {
  PREFIX=$(basename "$1" .fa)
  scripts/00_filter_helixer_annotation.py ./Helixer/00_raw/GFF/${PREFIX}_Helixer.gff \
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

The filtering criteria are as follows.

- Gene models containing an exon shorter than 3 bp are removed.
- `five_prime_UTR`, `three_prime_UTR` and `exon` features are excluded from the GFF, so that only
  gene, mRNA and CDS features are retained.
- For the secretome set, only models present in `${PREFIX}_Helixer_secretome.fa` are retained.

## 2.5. Output layout

```
Helixer/00_raw/GFF/                 # Raw Helixer GFF (PREFIX_Helixer.gff)
Helixer/00_raw/Proteomes/           # Raw translated proteins (gffread)
Helixer/00_raw/Transcriptomes/      # Raw CDS (gffread)
Helixer/00_raw/Secretomes/          # Raw Helixer-predicted secretome FASTA (input to the filtering script)
Helixer/10_gff_qc/                  # Quality-control-annotated Helixer GFF and logs
Helixer/20_filtered/GFF/            # Filtered Helixer proteome and secretome GFF
Helixer/20_filtered/Proteomes/      # Filtered Helixer proteome FASTA
Helixer/20_filtered/Secretomes/     # Filtered Helixer secretome FASTA
```
