# 3. Secretome prediction

**Objective:** identify the secreted proteins among the predicted proteomes of
[step 1](01_braker_prediction.md) and [step 2](02_helixer_prediction.md).

**Output:** `Secretomes/${PREFIX}_secretome.fa` and
`Helixer/00_raw/Secretomes/${PREFIX}_Helixer_secretome.fa`.

Secreted proteins were predicted from each proteome following the PexFinder procedure of
[Torto et al. (2003)](https://doi.org/10.1101/gr.910003): a signal peptide is required, and
proteins carrying a transmembrane domain or a mitochondrial targeting signal are removed.

```bash
# Signal peptide prediction (SignalP 2.0)
signalp -t euk -f summary -trunc 70 ${PREFIX}.proteins.fa

# Transmembrane domains in the proteins that passed (TMHMM 2.0)
tmhmm -short signalp_predicted_pex_proteins.fasta

# Mitochondrial targeting signals in the same proteins (TargetP 2.0)
targetp -batch 500 -fasta signalp_predicted_pex_proteins.fasta -prefix ${PREFIX} -tmp ${TMPDIR}
```

- The three predictors were run from a Singularity container.
- `-t euk`: eukaryotic SignalP model.
- `-trunc 70`: consider only the first 70 residues.
- A protein is retained when SignalP predicts a signal peptide with an HMM score of at least 0.8.
- A protein is discarded when TMHMM reports one or more predicted transmembrane helices
  (`PredHel`), or when TargetP assigns a mitochondrial targeting signal.

The proteins passing all three filters constitute the secretome of that proteome. The procedure was
applied to the BRAKER proteomes of [step 1](01_braker_prediction.md), giving
`Secretomes/${PREFIX}_secretome.fa`, and to the Helixer proteomes of
[step 2.3](02_helixer_prediction.md#23-extraction-of-protein-and-cds-sequences), giving
`Helixer/00_raw/Secretomes/${PREFIX}_Helixer_secretome.fa`. The latter is the input to the Helixer
filtering of [step 2.4](02_helixer_prediction.md#24-filtering-of-the-helixer-annotation), which
restricts it to the gene models retained there.

The BRAKER secretome and the filtered Helixer secretome are used to build the miniprot query set in
[step 4](04_reference_protein_datasets.md), and the BRAKER secretome is used again in
[step 11](11_helixer_specific_secreted_proteins.md).
