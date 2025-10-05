
### 8. Merge all evidence and rename IDs (60_merged)

**Objective:** combine BRAKER, Helixer secretome non‑overlap, DIAMOND‑retained Helixer proteome, miniprot rescues, then give stable IDs.

```bash
# Pack Helixer evidence into one temp GFF
gffread -M -K \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff \
  results/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff \
  > results/60_merged/temp/${PREFIX}_helixer_temp.gff

# Merge: BRAKER + Helixer (temp) + miniprot cleaned
gffread --sort-alpha --cluster-only --force-exons \
  results/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  results/60_merged/temp/${PREFIX}_helixer_temp.gff \
  results/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff \
  > results/60_merged/temp/${PREFIX}_temp.gff
```

- `--sort-alpha`: Sorts features alphabetically by sequence name and feature ID, ensuring a consistent order in the merged GFF.
- `--cluster-only`: Clusters transcripts into loci without merging or discarding any, preserving all input features.
- `--force-exons`: Ensures that exon features are explicitly included for all transcripts in the output, even if missing in the input.

```bash
# Rename to stable per‑project IDs
./70_rename_ids.py ${PREFIX} results/60_merged/temp/${PREFIX}_temp.gff \
  > results/60_merged/${PREFIX}_merged.gff
```

**About `70_rename_ids.py`:**

This script takes a merged GFF file and systematically renames all gene, transcript, exon, and CDS IDs to a consistent, project-specific format.  
- Each gene is assigned a unique identifier: `${PREFIX}_g00001`, `${PREFIX}_g00002`, etc.
- Transcripts are named as `${PREFIX}_g00001.t1`, `${PREFIX}_g00001.t2`, etc.
- Exons and CDS features are named as `${PREFIX}_g00001.t1.exon1`, `${PREFIX}_g00001.t1.CDS1`, etc.
