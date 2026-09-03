# 4. Reference protein datasets

**Objective:** assemble the protein sets that are used as queries in the miniprot searches
([step 5](05_miniprot_searches.md)) and in the DIAMOND search against characterized 70-15
proteins ([step 10](10_helixer_proteome_homology.md)).

**Output:** `all_secretomes.uniq.fa` and `data/70-15_RefSeq_characterized_protein.fa`.

## Contents

- [4.1. Secretome dataset of *M. oryzae* 70-15](#41-secretome-dataset-of-m-oryzae-70-15)
- [4.2. Published effector datasets](#42-published-effector-datasets)
- [4.3. Merging and deduplication of the query set](#43-merging-and-deduplication-of-the-query-set)
- [4.4. Characterized proteins of 70-15](#44-characterized-proteins-of-70-15)

## 4.1. Secretome dataset of *M. oryzae* 70-15

The secretome of *M. oryzae* strain 70-15 (assembly MG8) was extracted from the reference
proteome using the list of secreted proteins published by
[Yan et al. (2023)](https://doi.org/10.1093/plcell/koad036).

```bash
cat data/effector_datasets/Yan_et_al_2023_Supplemental_Data_Set_S7.txt | \
  samtools faidx -r - Magnaporthe_oryzae.MG8.pep.all.fa \
  > data/effector_datasets/Moryzae_MG8_XiaYan_secretome.fa
```

- `data/effector_datasets/Yan_et_al_2023_Supplemental_Data_Set_S7.txt` holds the 1,826 protein
  identifiers of the secretome.
- `samtools faidx -r`: extracts those identifiers from the MG8 proteome
  (`Magnaporthe_oryzae.MG8.pep.all.fa`, Ensembl Fungi).

## 4.2. Published effector datasets

The following experimentally characterized effector datasets were included in the query set.

| File | Sequences | Source |
|---|---|---|
| `data/effector_datasets/zif_effectors.fa` | 95 | ZiF effector family, [De la Concepcion et al., 2024](https://doi.org/10.1371/journal.ppat.1012277) |
| `data/effector_datasets/Bentham_et_al_2021_journal.ppat.1009957.s006.with_AVR-PikD_O23.fa` | 180 | APikL effector family, [Bentham et al., 2021](https://doi.org/10.1371/journal.ppat.1009957), together with the AVR-Pik_O23 sequence reported in [Sugihara et al., 2023](https://doi.org/10.1371/journal.pbio.3001945) |
| `data/effector_datasets/AVR-Mgk1_variants.fa` | 16 | AVR-Mgk1 reported in [Sugihara et al., 2023](https://doi.org/10.1371/journal.pbio.3001945) and its homologues |

## 4.3. Merging and deduplication of the query set

The 70-15 secretome, the published effector datasets and the predicted secretomes of all isolates
were concatenated and deduplicated.

```bash
cat Secretomes/*.fa Helixer/20_filtered/Secretomes/*.fa data/effector_datasets/*.fa \
  > all_secretomes.fa
seqkit rmdup -s -i -P all_secretomes.fa > all_secretomes.uniq.fa
```

- `Secretomes/*.fa`: secreted proteins predicted from the BRAKER gene models ([step 1](01_braker_prediction.md)).
- `Helixer/20_filtered/Secretomes/*.fa`: secreted proteins predicted from the filtered Helixer
  gene models ([step 2](02_helixer_prediction.md)).
- `data/effector_datasets/*.fa`: the 70-15 secretome of
  [4.1](#41-secretome-dataset-of-m-oryzae-70-15) and the effector datasets of
  [4.2](#42-published-effector-datasets).
- `-s`: remove duplicates on the basis of the sequence rather than the identifier.
- `-i`: ignore case.
- `-P`: consider only the positive strand when comparing sequences.

## 4.4. Characterized proteins of 70-15

Entries annotated as uncharacterized or hypothetical proteins were removed from the 70-15 RefSeq
proteome, so that only proteins with a functional annotation were used as queries.

```bash
scripts/01_remove_uncharacterized_seq.py data/70-15_RefSeq_protein.fa \
  > data/70-15_RefSeq_characterized_protein.fa
```

`01_remove_uncharacterized_seq.py` discards any record whose description contains
`uncharacterized protein` or `hypothetical protein` and writes the remaining records to standard
output.
