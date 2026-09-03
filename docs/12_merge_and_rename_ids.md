# 12. Merging of all evidence and renaming of identifiers

**Objective:** combine the BRAKER models with the non-redundant Helixer and miniprot models
selected in [step 9](09_remove_overlapping_models.md), [step 10](10_helixer_proteome_homology.md)
and [step 11](11_helixer_specific_secreted_proteins.md), and assign systematic identifiers to the
merged annotation.

**Output:** `${OUTDIR}/60_merged/${PREFIX}_merged.gff`.

```bash
# Combine the three Helixer-derived sets into a single temporary GFF
gffread -M -K \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff \
  ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff \
  > ${OUTDIR}/60_merged/temp/${PREFIX}_helixer_temp.gff

# Merge BRAKER, Helixer and miniprot models
gffread --sort-alpha --cluster-only --force-exons \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
  ${OUTDIR}/60_merged/temp/${PREFIX}_helixer_temp.gff \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff \
  > ${OUTDIR}/60_merged/temp/${PREFIX}_temp.gff
```

- `-M -K`: cluster the three Helixer-derived sets into loci and discard redundant transcripts,
  which can arise because the same model may be selected in more than one of steps 9, 10 and 11.
- `--sort-alpha`: sort the features alphabetically by sequence name and feature identifier, so
  that the merged file has a deterministic order.
- `--cluster-only`: cluster the transcripts into loci without merging or discarding any of them,
  so that all models selected in the preceding steps are preserved.
- `--force-exons`: write explicit exon features for all transcripts, including those for which the
  input contains only CDS features.

```bash
# Assign systematic identifiers
scripts/70_rename_ids.py ${PREFIX} ${OUTDIR}/60_merged/temp/${PREFIX}_temp.gff \
  > ${OUTDIR}/60_merged/${PREFIX}_merged.gff
```

## `70_rename_ids.py`

The script renames all gene, transcript, exon and CDS identifiers of the merged GFF file to a
consistent, isolate-specific format.

| Feature | Format | Example |
|---|---|---|
| Gene | `${PREFIX}_gNNNNN` | `${PREFIX}_g00001` |
| Transcript | `${PREFIX}_gNNNNN.tN` | `${PREFIX}_g00001.t1` |
| Exon | `${PREFIX}_gNNNNN.tN.exonN` | `${PREFIX}_g00001.t1.exon1` |
| CDS | `${PREFIX}_gNNNNN.tN.CDSN` | `${PREFIX}_g00001.t1.CDS1` |

Genes are numbered consecutively in the order in which they appear in the sorted, merged file.
