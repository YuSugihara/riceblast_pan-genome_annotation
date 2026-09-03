# Command reference

This document collects the options of the third-party tools used in the pipeline, with the step in
which each is used.

## Extraction of features from a GFF3 file

Extract all mRNA features ([step 9](09_remove_overlapping_models.md)):

```bash
grep 'mRNA' input.gff > output.mrna.gff
```

Extract the `ID` attribute from the ninth column ([step 9](09_remove_overlapping_models.md)):

```bash
cut -f9 file.gff | cut -d ';' -f1 | cut -d '=' -f2 > ids.txt
```

Extract the features listed in a file of identifiers ([steps 9](09_remove_overlapping_models.md),
[10](10_helixer_proteome_homology.md) and [11](11_helixer_specific_secreted_proteins.md)):

```bash
gffread --ids ids.txt input.gff > output.filtered.gff
```

- `--ids`: file containing one identifier per line.

Extract features while clustering the transcripts into loci and discarding redundant transcripts
([steps 9](09_remove_overlapping_models.md) and [12](12_merge_and_rename_ids.md)):

```bash
gffread -M -K --ids ids.txt input.gff > output.filtered.gff
```

- `-M`, `--merge`: cluster the input transcripts into loci and discard redundant transcripts,
  that is, transcripts with identical intron chains whose boundaries are equal to, or contained
  within, those of another transcript.
- `-K`: in combination with `-M`, also discard shorter transcripts that are fully contained within
  another transcript whose intron chain matches theirs, so that the most representative transcript
  of each locus is retained.

## Removal of overlapping features

Remove any feature of file A that overlaps a feature of file B
([step 9](09_remove_overlapping_models.md)):

```bash
bedtools subtract -A -a fileA.gff -b fileB.gff > fileA.nonoverlap.gff
```

- `-A`: remove the overlapping feature of A in its entirety, rather than only the overlapping
  interval.

## Extraction of sequences from a GFF3 file and a genome

Extract CDS sequences ([steps 6](06_extract_cds.md) and [13](13_export_final_sequences.md)):

```bash
gffread -x output.cds.fa -g genome.fa input.gff
```

Extract protein sequences ([steps 2](02_helixer_prediction.md),
[10](10_helixer_proteome_homology.md), [11](11_helixer_specific_secreted_proteins.md) and
[13](13_export_final_sequences.md)):

```bash
gffread -y output.protein.fa -g genome.fa input.gff
```

- `-x`: write CDS sequences.
- `-y`: write protein sequences.
- `-g`: reference genome from which the sequences are extracted.

## Deduplication of a FASTA file

Remove duplicate sequences ([step 4](04_reference_protein_datasets.md)):

```bash
seqkit rmdup -s -i -P input.fa > output.uniq.fa
```

- `-s`: remove duplicates on the basis of the sequence rather than the identifier.
- `-i`: ignore case.
- `-P`: consider only the positive strand.

## Merging and sorting of GFF3 files

Merge, sort and cluster GFF3 files while writing explicit exon features
([step 12](12_merge_and_rename_ids.md)):

```bash
gffread --sort-alpha --cluster-only --force-exons input1.gff input2.gff > merged.gff
```

- `--sort-alpha`: sort the features alphabetically by sequence name and feature identifier.
- `--cluster-only`: cluster the transcripts into loci without merging or discarding any of them.
- `--force-exons`: write explicit exon features for all transcripts.

## Auxiliary scripts

The following scripts are not invoked by `scripts/annotation_pipeline_script.sh`. They were used to
inspect the intermediate results.

### `60_get_same_seq_id.py`

```bash
scripts/60_get_same_seq_id.py query.fa subject.fa
```

Reports the identifiers of the records of `subject.fa` whose name begins with `SEC`, that is,
miniprot models derived from the secretome query panel ([step 5](05_miniprot_searches.md)), and
whose protein sequence is identical to a sequence present in `query.fa`. It was used to check
whether miniprot-derived secreted proteins are already represented in another protein set.

### `61_check_effector_seq.py`

```bash
scripts/61_check_effector_seq.py gffcompare.loci braker_proteome.fa helixer_secretome.fa miniprot.fa
```

Reads a gffcompare `.loci` file together with the BRAKER proteome, the Helixer secretome and the
miniprot protein sequences, and reports miniprot transcripts that share a locus with a BRAKER or
Helixer model but whose protein sequence is not identical to any of them. It was used to inspect
the residual redundancy between the miniprot models and the other evidence tracks.
