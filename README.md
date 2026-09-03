# Rice blast pan-genome annotation pipeline

<img src="figures/Figure_1_Annotation_pipeline.png" width=800>

**Figure 1.** Overview of the annotation pipeline.

This repository documents the annotation workflow used to construct a non-redundant gene set for
each *Magnaporthe oryzae* genome assembly analysed in this study. It is intended as a record of
the procedure that accompanies the Materials and Methods of the associated publication, and it
provides the exact commands, parameters and helper scripts that were used, together with the
input and output layout of each step.

The workflow integrates three sources of evidence: gene prediction with BRAKER, trained on
RNA-seq alignments and protein hints; deep-learning-based prediction with Helixer; and
protein-to-genome alignment with miniprot.
Models from these sources are quality-filtered, resolved against each other by a fixed source
priority, merged, and renamed, yielding a GFF3 file together with CDS and protein FASTA files
for each assembly.

## Contents

| Step | Document | Output directory |
|---|---|---|
| – | [Software dependencies](docs/00_dependencies.md) | – |
| 1 | [Gene prediction with BRAKER](docs/01_braker_prediction.md) | `Gene_Models/`, `Proteomes/` |
| 2 | [Gene prediction with Helixer](docs/02_helixer_prediction.md) | `Helixer/` |
| 3 | [Secretome prediction](docs/03_secretome_prediction.md) | `Secretomes/`, `Helixer/00_raw/Secretomes/` |
| 4 | [Reference protein datasets](docs/04_reference_protein_datasets.md) | `all_secretomes.uniq.fa` |
| 5 | [miniprot searches](docs/05_miniprot_searches.md) | `results/00_miniprot/` |
| 6 | [CDS extraction from each source](docs/06_extract_cds.md) | `results/10_cds/` |
| 7 | [Quality control of each GFF](docs/07_gff_qc.md) | `results/20_gff_qc/` |
| 8 | [Filtering of gene models](docs/08_filter_gene_models.md) | `results/30_filtered_gff/` |
| 9 | [Removal of overlapping gene models](docs/09_remove_overlapping_models.md) | `results/40_non_overlap_gff/` |
| 10 | [Homology screening of Helixer proteome models](docs/10_helixer_proteome_homology.md) | `results/40_non_overlap_gff/` |
| 11 | [Helixer-specific secreted proteins](docs/11_helixer_specific_secreted_proteins.md) | `results/50_helixer_uniq_effectors/` |
| 12 | [Merging of all evidence and renaming of IDs](docs/12_merge_and_rename_ids.md) | `results/60_merged/` |
| 13 | [Export of final sequences](docs/13_export_final_sequences.md) | `results/70_extracted_seq/` |
| – | [Command reference](docs/14_command_reference.md) | – |
| – | [Effector annotation benchmark](benchmark/README.md) | – |

The numbering of the documents follows the order of execution, and is independent of the numeric
prefixes of the output directories and of the helper scripts.

`scripts/annotation_pipeline_script.sh` builds the query sets of step 4 and then runs steps 5 to
13 for all assemblies in parallel. The individual documents describe the same commands one step at
a time.

The scripts and input datasets are grouped here into `scripts/` and `data/`, and the commands in
the documentation refer to these locations. In the working directory used for the analysis they
resided at the top level, as the paths inside `scripts/annotation_pipeline_script.sh` show. A
complete archive of that directory is available as Supplementary Data ####.

## Conventions

The following variables and directory names are used throughout the documentation.

| Name | Meaning |
|---|---|
| `${PREFIX}` | Isolate identifier, derived from the assembly file name (`basename ${PREFIX}.fa .fa`) |
| `${OUTDIR}` | Root of the pipeline output, `./results` |
| `Frozen_assemblies/${PREFIX}.fa` | Soft-masked genome assembly used as input for all steps |
| `Gene_Models/${PREFIX}.gff3` | BRAKER gene models (step 1) |
| `Proteomes/${PREFIX}.proteins.fa` | Protein sequences of the BRAKER gene models |
| `Secretomes/${PREFIX}_secretome.fa` | BRAKER-predicted secreted proteins |
| `Helixer/` | Raw, quality-controlled and filtered Helixer predictions (step 2) |
| `data/effector_datasets/` | Published effector and secretome datasets used as miniprot queries (step 4) |

## Input datasets

The files below are inputs that the pipeline did not produce. They are collected in `data/`.

| File | Description |
|---|---|
| `data/moryzae_repeats.fa.classified` | Custom library of known *M. oryzae* repeat elements, used for soft-masking in [step 1](docs/01_braker_prediction.md) |
| `data/effector_datasets/` | Published effector and secretome datasets used to build the miniprot query set in [step 4](docs/04_reference_protein_datasets.md) |
| `data/prefer_braker2.cfg` | TSEBRA configuration used to reconcile the two BRAKER predictions in [step 1](docs/01_braker_prediction.md) |
| `data/Mo_Effectors.fasta` | Curated set of known *M. oryzae* effectors and their homologues; one of the two components of the protein query set used for BRAKER in [step 1](docs/01_braker_prediction.md), the other being the 70-15 reference proteome (UniProt proteome UP000009058) |
| `data/70-15_RefSeq_protein.fa` | Complete RefSeq proteome of *M. oryzae* strain 70-15 |
| `data/70-15_RefSeq_characterized_protein.fa` | Subset of the above from which entries annotated as uncharacterized or hypothetical proteins have been removed by `01_remove_uncharacterized_seq.py`; used as a miniprot query set and as a DIAMOND query set |

## Scripts

| Script | Used in | Purpose |
|---|---|---|
| `scripts/00_filter_helixer_annotation.py` | [Step 2](docs/02_helixer_prediction.md) | Annotates raw Helixer GFF with exon-length and secretion status |
| `scripts/01_remove_uncharacterized_seq.py` | [Step 4](docs/04_reference_protein_datasets.md) | Removes uncharacterized and hypothetical proteins from a FASTA file |
| `scripts/20_gff_qc.py` | [Step 7](docs/07_gff_qc.md) | Annotates gene models with quality-control tags |
| `scripts/40_clean_up_miniprot_annotation.py` | [Step 9](docs/09_remove_overlapping_models.md) | Selects miniprot transcripts terminating in a valid stop codon |
| `scripts/41_filter_diamond_results.py` | [Step 10](docs/10_helixer_proteome_homology.md) | Filters DIAMOND hits by identity and query coverage |
| `scripts/50_check_overlapping_effectors.py` | [Step 11](docs/11_helixer_specific_secreted_proteins.md) | Identifies Helixer-specific secreted proteins at BRAKER loci |
| `scripts/70_rename_ids.py` | [Step 12](docs/12_merge_and_rename_ids.md) | Assigns systematic gene, transcript, exon and CDS identifiers |
| `scripts/60_get_same_seq_id.py`, `scripts/61_check_effector_seq.py` | – | Auxiliary scripts used for inspection of the results; see the [command reference](docs/14_command_reference.md) |
