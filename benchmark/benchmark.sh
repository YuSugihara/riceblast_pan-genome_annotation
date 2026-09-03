#!/usr/bin/env bash
set -euo pipefail

# Benchmark how completely the final merged annotations describe known effector genes.
#
# Annotation evidence is taken only from the user-provided merged GFF directory
# and the matching protein FASTA directory. No intermediate annotation files are
# used.
#
# One BLASTP search is run per isolate against the complete merged proteome. The
# cumulative source-based annotation sets (AUGUSTUS, AUGUSTUS+Helixer,
# AUGUSTUS+Helixer+miniprot) are derived in 01_build_locus_tables.py from the
# second column of the merged GFF. `AUGUSTUS` is the source value that BRAKER
# assigns to its gene models.
#
# All analysis thresholds are options of the Python scripts, not of this
# wrapper. Pass them after `--`; see --help.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ASSEMBLY_DIR=""
MERGED_GFF_DIR=""
ANNOTATION_PROTEIN_DIR=""
EFFECTOR_FASTA="${SCRIPT_DIR}/data/benchmark_effectors.fa"
OUTDIR="./results"
THREADS=4
JOBS=4
TBLASTN_EVALUE="1e-5"
TBLASTN_MAX_TARGET_SEQS=20
BLASTP_EVALUE="1e-5"
DIAMOND_MAX_TARGET_SEQS=25
PYTHON_CMD="python3"
TABLE_ARGS=()

usage() {
  cat <<'USAGE'
Usage: benchmark.sh --assembly-dir DIR --merged-gff-dir DIR --protein-dir DIR [options] [-- TABLE_OPTIONS]

Required:
  --assembly-dir DIR        genome assemblies named <isolate>.fa
  --merged-gff-dir DIR      merged annotations named <isolate>_merged.gff
  --protein-dir DIR         proteins of those annotations named <isolate>_protein.fa

Options:
  --effector-fasta FILE     benchmark effector queries (default: data/benchmark_effectors.fa)
  --outdir DIR              output directory (default: ./results)
  --threads N               threads per search job (default: 4)
  --jobs N                  isolates processed in parallel (default: 4)
                            peak CPU use is roughly --jobs x --threads
  --tblastn-evalue X        tblastn E-value cutoff (default: 1e-5)
  --tblastn-max-target-seqs N   tblastn -max_target_seqs (default: 20)
  --blastp-evalue X         DIAMOND blastp E-value cutoff (default: 1e-5)
  --diamond-max-target-seqs N   DIAMOND --max-target-seqs (default: 25)
  --python PATH             Python interpreter with matplotlib (default: python3)
  -h, --help                show this message

Everything after `--` is passed verbatim to scripts/01_build_locus_tables.py,
which holds the analysis thresholds. For example:

  benchmark.sh --assembly-dir A --merged-gff-dir G --protein-dir P \
    -- --model-qcov-min 90 --overextended-len-ratio 1.5

Run `scripts/01_build_locus_tables.py --help` for the full list.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly-dir) ASSEMBLY_DIR="$2"; shift 2 ;;
    --merged-gff-dir) MERGED_GFF_DIR="$2"; shift 2 ;;
    --protein-dir) ANNOTATION_PROTEIN_DIR="$2"; shift 2 ;;
    --effector-fasta) EFFECTOR_FASTA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --tblastn-evalue) TBLASTN_EVALUE="$2"; shift 2 ;;
    --tblastn-max-target-seqs) TBLASTN_MAX_TARGET_SEQS="$2"; shift 2 ;;
    --blastp-evalue) BLASTP_EVALUE="$2"; shift 2 ;;
    --diamond-max-target-seqs) DIAMOND_MAX_TARGET_SEQS="$2"; shift 2 ;;
    --python) PYTHON_CMD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; TABLE_ARGS=("$@"); break ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

QUERY_DIR="${OUTDIR}/00_query"
GENOME_DIR="${OUTDIR}/10_tblastn_genome"
PROTEIN_DIR="${OUTDIR}/20_annotation_proteins"
BLASTP_DIR="${OUTDIR}/30_blastp"
SUMMARY_DIR="${OUTDIR}/40_summary"
PLOT_DIR="${OUTDIR}/50_plots"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $1" >&2
    exit 1
  fi
}

require_dir() {
  if [[ -z "$2" ]]; then
    echo "ERROR: $1 is required" >&2
    usage >&2
    exit 2
  fi
  if [[ ! -d "$2" ]]; then
    echo "ERROR: directory not found: $2" >&2
    exit 1
  fi
}

check_inputs() {
  require_command seqkit
  require_command diamond
  require_command makeblastdb
  require_command tblastn

  require_dir --assembly-dir "${ASSEMBLY_DIR}"
  require_dir --merged-gff-dir "${MERGED_GFF_DIR}"
  require_dir --protein-dir "${ANNOTATION_PROTEIN_DIR}"

  if [[ ! -s "${EFFECTOR_FASTA}" ]]; then
    echo "ERROR: benchmark effector FASTA not found: ${EFFECTOR_FASTA}" >&2
    exit 1
  fi
}

list_genomes() {
  find "${ASSEMBLY_DIR}" -maxdepth 1 -type f -name '*.fa' ! -name '._*' | sort
}

check_isolates() {
  local missing=0 genome prefix
  while read -r genome; do
    prefix="$(basename "${genome}" .fa)"
    if [[ ! -s "${MERGED_GFF_DIR}/${prefix}_merged.gff" ]]; then
      echo "ERROR: missing merged GFF for ${prefix}: ${MERGED_GFF_DIR}/${prefix}_merged.gff" >&2
      missing=1
    fi
    if [[ ! -s "${ANNOTATION_PROTEIN_DIR}/${prefix}_protein.fa" ]]; then
      echo "ERROR: missing protein FASTA for ${prefix}: ${ANNOTATION_PROTEIN_DIR}/${prefix}_protein.fa" >&2
      missing=1
    fi
  done < <(list_genomes)
  [[ "${missing}" -eq 0 ]] || exit 1
}

prepare_query() {
  mkdir -p "${QUERY_DIR}"
  cp "${EFFECTOR_FASTA}" "${QUERY_DIR}/benchmark_effectors.fa"
  seqkit fx2tab -n -l "${QUERY_DIR}/benchmark_effectors.fa" > "${QUERY_DIR}/benchmark_effectors.length.tsv"
}

prepare_annotation_proteins() {
  # gffread writes a trailing '.' for the stop codon and '.' for unresolved
  # residues; DIAMOND rejects '.' in a protein sequence, so replace it with 'X'.
  # Sequence identifiers are left untouched so that they match the mRNA IDs of
  # the merged GFF.
  local prefix="$1"
  awk '
    /^>/ { print; next }
    { gsub(/\./, "X"); print }
  ' "${ANNOTATION_PROTEIN_DIR}/${prefix}_protein.fa" \
    > "${PROTEIN_DIR}/${prefix}/${prefix}.protein.fa"
}

run_genome_tblastn() {
  local genome="$1"
  local prefix="$2"
  local db="${GENOME_DIR}/${prefix}/${prefix}.genome.db"

  makeblastdb -in "${genome}" -dbtype nucl -out "${db}" \
    > "${GENOME_DIR}/${prefix}/${prefix}.makeblastdb.log" 2>&1
  tblastn \
    -query "${QUERY_DIR}/benchmark_effectors.fa" \
    -db "${db}" \
    -evalue "${TBLASTN_EVALUE}" \
    -max_target_seqs "${TBLASTN_MAX_TARGET_SEQS}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" \
    -num_threads "${THREADS}" \
    -out "${GENOME_DIR}/${prefix}/${prefix}.tblastn.tsv" \
    2> "${GENOME_DIR}/${prefix}/${prefix}.tblastn.log"
}

run_blastp() {
  local prefix="$1"
  local protein_fasta="${PROTEIN_DIR}/${prefix}/${prefix}.protein.fa"
  local db="${BLASTP_DIR}/${prefix}/${prefix}.db"

  diamond makedb --threads 1 --in "${protein_fasta}" -d "${db}" \
    > "${BLASTP_DIR}/${prefix}/${prefix}.makedb.log" 2>&1
  diamond blastp \
    --threads "${THREADS}" \
    --very-sensitive \
    --evalue "${BLASTP_EVALUE}" \
    --max-target-seqs "${DIAMOND_MAX_TARGET_SEQS}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp qlen slen \
    -d "${db}" \
    -q "${QUERY_DIR}/benchmark_effectors.fa" \
    -o "${BLASTP_DIR}/${prefix}/${prefix}.blastp.tsv" \
    2> "${BLASTP_DIR}/${prefix}/${prefix}.blastp.log"
}

run_one_isolate() {
  local genome="$1"
  local prefix
  prefix="$(basename "${genome}" .fa)"
  echo "${prefix}"
  mkdir -p "${GENOME_DIR}/${prefix}" "${PROTEIN_DIR}/${prefix}" "${BLASTP_DIR}/${prefix}"

  run_genome_tblastn "${genome}" "${prefix}"
  prepare_annotation_proteins "${prefix}"
  run_blastp "${prefix}"
}

build_tables() {
  "${PYTHON_CMD}" "${SCRIPT_DIR}/scripts/01_build_locus_tables.py" \
    --outdir "${OUTDIR}" \
    --assembly-dir "${ASSEMBLY_DIR}" \
    --merged-gff-dir "${MERGED_GFF_DIR}" \
    ${TABLE_ARGS[@]+"${TABLE_ARGS[@]}"}
}

plot_results() {
  "${PYTHON_CMD}" "${SCRIPT_DIR}/scripts/02_plot_effector_benchmark.py" \
    --summary-dir "${SUMMARY_DIR}" \
    --plot-dir "${PLOT_DIR}"
}

main() {
  check_inputs

  local genome_count
  genome_count="$(list_genomes | wc -l | tr -d ' ')"
  if [[ "${genome_count}" -eq 0 ]]; then
    echo "ERROR: no *.fa assemblies found in ${ASSEMBLY_DIR}" >&2
    exit 1
  fi
  echo "Found ${genome_count} assemblies in ${ASSEMBLY_DIR}"
  check_isolates

  prepare_query

  export ANNOTATION_PROTEIN_DIR QUERY_DIR GENOME_DIR PROTEIN_DIR BLASTP_DIR
  export THREADS TBLASTN_EVALUE TBLASTN_MAX_TARGET_SEQS
  export BLASTP_EVALUE DIAMOND_MAX_TARGET_SEQS
  export -f prepare_annotation_proteins run_genome_tblastn run_blastp run_one_isolate

  list_genomes \
    | xargs -P "${JOBS}" -I{} bash -c 'run_one_isolate "$@"' _ {}

  build_tables
  plot_results
}

main
