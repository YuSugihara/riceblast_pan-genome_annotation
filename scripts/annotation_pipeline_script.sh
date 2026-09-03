cat ./effector_datasets/Yan_et_al_2023_dataset/Yan_et_al_2023_Supplemental_Data_Set_S7.txt | \
samtools faidx -r - \
               ./effector_datasets/Yan_et_al_2023_dataset/Magnaporthe_oryzae.MG8.pep.all.fa > \
               ./effector_datasets/Moryzae_MG8_XiaYan_secretome.fa


cat Secretomes/*.fa \
    Helixer/20_filtered/Secretomes/*.fa \
    effector_datasets/*.fa > \
    all_secretomes.fa

seqkit rmdup -s \
             -i \
             -P \
             all_secretomes.fa > \
             all_secretomes.uniq.fa

./01_remove_uncharacterized_seq.py 70-15_RefSeq_protein.fa > \
                                   70-15_RefSeq_characterized_protein.fa

function run_annotation_pipeline () {
  PREFIX=`basename $1 .fa`
  OUTDIR=./results
  echo ${PREFIX}

  mkdir -p ${OUTDIR}/00_miniprot/${PREFIX}
  mkdir -p ${OUTDIR}/10_cds/${PREFIX}
  mkdir -p ${OUTDIR}/20_gff_qc/${PREFIX}
  mkdir -p ${OUTDIR}/30_filtered_gff/${PREFIX}
  mkdir -p ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond
  mkdir -p ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}
  mkdir -p ${OUTDIR}/60_merged/temp
  mkdir -p ${OUTDIR}/70_extracted_seq/protein
  mkdir -p ${OUTDIR}/70_extracted_seq/cds

  ## 00_miniprot: Run miniprot

  miniprot -t 2 \
           -G 3k \
           --gff \
           -P SEC \
           --outs=0.5 \
           -p 0.3 \
           $1 \
           all_secretomes.uniq.fa \
           2> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "##PAF" > ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff

  miniprot -t 2 \
           -G 3k \
           --gff \
           -P RefSeq \
           $1 \
           70-15_RefSeq_characterized_protein.fa \
           2>> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.log | \
  grep -v "#" >> ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff

  ## 10_cds: Extract CDS sequences

  gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.braker.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ./Gene_Models/${PREFIX}.gff3

  gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff

  gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff

  gffread -x ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff

  ## 20_gff_qc: GFF3 QC

  ./20_gff_qc.py ./Gene_Models/${PREFIX}.gff3 \
                 ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.braker.fa \
                 150,180,195 \
                 10,25,50 \
                 1> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.braker_qc.gff \
                 2> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.braker_qc.txt

  ./20_gff_qc.py ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.proteome_filtered.gff \
                 ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa \
                 150,180,195 \
                 10,25,50 \
                 1> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_proteome_qc.gff \
                 2> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_proteome_qc.txt

  ./20_gff_qc.py ./Helixer/20_filtered/GFF/${PREFIX}.helixer_pre_qc.secretome_filtered.gff \
                 ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa \
                 150,180,195 \
                 10,25,50 \
                 1> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_secretome_qc.gff \
                 2> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_secretome_qc.txt

  ./20_gff_qc.py ${OUTDIR}/00_miniprot/${PREFIX}/${PREFIX}.miniprot.gff \
                 ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa \
                 150,180,195 \
                 10,25,50 \
                 1> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.miniprot_qc.gff \
                 2> ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.miniprot_qc.txt

  ## 30_filtered_gff: Filter GFF3

  grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.braker_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff

  grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_proteome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff

  grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.helixer_secretome_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'no_stop_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff

  grep -v 'not_multiple_of_3' ${OUTDIR}/20_gff_qc/${PREFIX}/${PREFIX}.miniprot_qc.gff | \
  grep -v 'stop_codon_in_cds' | \
  grep -v 'no_start_codon' | \
  grep -v 'shorter_than_150nt' | \
  grep -v 'masked_over_25' | \
  cut -f 1-9 > \
  ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff

  ## 40_non_overlap_gff: Remove overlapping genes
  
  ## Remove overlapping genes for Helixer secretome
  grep 'mRNA' \
       ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff > \
       ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff

  bedtools subtract -A \
                    -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff \
                    -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff > \
                    ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff

  cut -f 9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f 1 | \
  cut -d '=' -f 2 > \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt

  ## Remove overlapping genes for miniprot

  grep 'mRNA' \
       ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff > \
       ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff

  bedtools subtract -A \
                    -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff \
                    -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
  bedtools subtract -A \
                    -a stdin \
                    -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff > \
                    ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff

  cut -f 9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f 1 | \
  cut -d '=' -f 2 > \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt

  ## Remove overlapping genes for Helixer proteome

  grep 'mRNA' \
       ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff > \
       ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff

  bedtools subtract -A \
                    -a ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff \
                    -b ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff | \
  bedtools subtract -A \
                    -a stdin \
                    -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff | \
  bedtools subtract -A \
                    -a stdin \
                    -b ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff > \
                    ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff

  cut -f 9 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff | \
  cut -d ';' -f 1 | \
  cut -d '=' -f 2 > \
  ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt

  ## Extract non-overlapping GFF lines

  gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.txt \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff > \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff

  gffread -M \
          -K \
          --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.txt \
                ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff > \
                ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff

  ./40_clean_up_miniprot_annotation.py ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.gff \
                                       ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa > \
                                       ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt

  gffread -M \
          -K \
          --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.txt \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.gff > \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff

  gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.txt \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff > \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff


  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.gff
  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.mRNA.non_overlap.gff
  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.gff
  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.mRNA.non_overlap.gff
  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.gff
  rm -f ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.mRNA.non_overlap.gff

  ## Check conserved proteins annotated by Helixer

  gffread -y ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.gff

  diamond makedb --threads 1 \
                 --in ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
                 -d ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
                 2> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.log

  diamond blastp --threads 1 \
                 --very-sensitive \
                 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
                 -d ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
                 -q 70-15_RefSeq_characterized_protein.fa \
                 -o ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv \
                 2>> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.log

  ./41_filter_diamond_results.py ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa.tsv \
                                 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.faa \
                                 1> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv  \
                                 2> ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa

  ## Extract Helixer proteins that are related to characterized proteins of 70-15

  cut -f 2 ${OUTDIR}/40_non_overlap_gff/${PREFIX}/diamond/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.filtered.faa.tsv | \
  sort | \
  uniq > ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt

  gffread --ids ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.txt \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.gff > \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff

  # 50_helixer_uniq_effectors: Extract Helixer unique effectors which overlap with BRAKER annotation

  gffcompare -p ${PREFIX} \
             -T \
             -o ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare \
             ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
             ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff \
             2> /dev/null

  ./50_check_overlapping_effectors.py ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_50_gffcompare.loci \
                                      ./Proteomes/${PREFIX}.proteins.fa \
                                      ./Secretomes/${PREFIX}_secretome.fa \
                                      ./Helixer/20_filtered/Proteomes/${PREFIX}_Helixer.proteome_filtered.fa \
                                      ./Helixer/20_filtered/Secretomes/${PREFIX}_Helixer.secretome_filtered.fa \
                                      1> ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.tsv \
                                      2> ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt

  gffread --ids ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.txt \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.gff > \
          ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff

  gffread -y ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff

    ## 60_merged: Merge GFF3

  gffread -M \
          -K \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_secretome_qc.filtered.non_overlap.gff \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.helixer_proteome_qc.filtered.non_overlap.70-15_refseq.gff \
          ${OUTDIR}/50_helixer_uniq_effectors/${PREFIX}/${PREFIX}_helixer_uniq_effectors.gff > \
          ${OUTDIR}/60_merged/temp/${PREFIX}_helixer_temp.gff

  gffread --sort-alpha \
          --cluster-only \
          --force-exons \
          ${OUTDIR}/30_filtered_gff/${PREFIX}/${PREFIX}.braker_qc.filtered.gff \
          ${OUTDIR}/60_merged/temp/${PREFIX}_helixer_temp.gff \
          ${OUTDIR}/40_non_overlap_gff/${PREFIX}/${PREFIX}.miniprot_qc.filtered.non_overlap.cleaned_up.gff > \
          ${OUTDIR}/60_merged/temp/${PREFIX}_temp.gff

  ./70_rename_ids.py ${PREFIX} \
                     ${OUTDIR}/60_merged/temp/${PREFIX}_temp.gff > \
                     ${OUTDIR}/60_merged/${PREFIX}_merged.gff


  # 70_extracted_seq: Extract protein and CDS sequences

  gffread -y ${OUTDIR}/70_extracted_seq/protein/${PREFIX}_protein.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ${OUTDIR}/60_merged/${PREFIX}_merged.gff

  gffread -x ${OUTDIR}/70_extracted_seq/cds/${PREFIX}_cds.fa \
          -g ./Frozen_assemblies/${PREFIX}.fa \
          ${OUTDIR}/60_merged/${PREFIX}_merged.gff

  gzip ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.braker.fa
  gzip ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_proteome.fa
  gzip ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.helixer_secretome.fa
  gzip ${OUTDIR}/10_cds/${PREFIX}/${PREFIX}.miniprot.fa

}

export -f run_annotation_pipeline

find ./Frozen_assemblies/*.fa | \
xargs -P 4 -I"%" bash -c "run_annotation_pipeline %"

