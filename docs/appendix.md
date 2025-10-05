## Appendix

### 1. GFF and annotation file manipulation tips

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

#### Merging and sorting GFF files with clustering

- **Merge, sort, and cluster GFF files, ensuring all exons are present:**
  ```bash
  gffread --sort-alpha --cluster-only --force-exons input1.gff input2.gff > merged.gff
  ```
  - `--sort-alpha`: Sorts features alphabetically by sequence name and feature ID.
  - `--cluster-only`: Clusters transcripts into loci without merging or discarding any.
  - `--force-exons`: Ensures exon features are explicitly included for all transcripts.
