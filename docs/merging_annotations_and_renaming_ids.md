## Merging annotations and renaming IDs


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
| Helixer           | 0.3.2    | [Helixer](https://github.com/weberlab-hhu/Helixer) |
| RepeatMasker     | 4.0.9    | [RepeatMasker](http://www.repeatmasker.org/) |

---

## Table of Contents

  - [0. Build secretome dataset](docs/0_build_secretome_dataset.md)
  - [1. miniprot searches (00_miniprot)](docs/1_miniprot_searches.md)
  - [2. Extract CDS from each source (10_cds)](docs/2_extract_cds_from_each_source.md)
  - [3. QC each GFF (20_gff_qc)](docs/3_qc_each_gff.md)
  - [4. Filter gene models (30_filtered_gff)](docs/4_filter_gene_models.md)
  - [5. Remove overlaps to define non-redundant gene sets (40_non_overlap_gff)](docs/5_remove_overlaps_to_define_non_redundant_gene_sets.md)
  - [6. Retain characterized proteins (Helixer ↔ 70-15 RefSeq)](docs/6_retain_characterized_proteins.md)
  - [7. Helixer-unique secreted proteins overlapping BRAKER (50_helixer_uniq_secreted_proteins)](docs/7_helixer_uniq_secreted_proteins.md)
  - [8. Merge all evidence and rename IDs (60_merged)](docs/8_merge_all_evidence_and_rename_ids.md)
  - [9. Export final sequences (70_extracted_seq)](docs/9_export_final_sequences.md)

