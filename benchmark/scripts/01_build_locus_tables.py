#!/usr/bin/env python3
"""Build locus-level benchmark tables for effector annotation.

Annotation evidence comes only from the merged GFF files, which supply the gene
coordinates and the GFF source values, and from the matching protein FASTA files
that the wrapper searches with BLASTP.

The cumulative annotation sets are defined from the GFF source column:
AUGUSTUS, AUGUSTUS+Helixer, and AUGUSTUS+Helixer+miniprot. `AUGUSTUS` is the
source value that BRAKER assigns to its gene models. One BLASTP search is run
per isolate against the whole merged proteome, and hits are assigned to the sets
here, from the GFF source of the transcript each hit belongs to.

All thresholds are command-line options.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}
STAGES = ("AUGUSTUS", "AUGUSTUS_Helixer", "Final")
STAGE_SOURCES = {
    "AUGUSTUS": {"AUGUSTUS"},
    "AUGUSTUS_Helixer": {"AUGUSTUS", "Helixer"},
    "Final": {"AUGUSTUS", "Helixer", "miniprot"},
}


@dataclass
class Config:
    outdir: Path
    assembly_dir: Path
    merged_gff_dir: Path
    genome_qcov_min: float
    genome_complete_qcov_min: float
    hard_genome_qcov_min: float
    effector_pident_min: float
    hard_effector_pident_min: float
    genomic_orf_aa_min: int
    locus_merge_distance: int
    model_qcov_min: float
    model_len_ratio_min: float
    overextended_len_ratio: float


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate_nt(seq: str) -> str:
    seq = seq.upper()
    return "".join(GENETIC_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def aa_identity(query: str, subject: str) -> float:
    if not query or not subject:
        return 0.0
    n = min(len(query), len(subject))
    return 100.0 * sum(query[i] == subject[i] for i in range(n)) / len(query)


def parse_attrs(value: str) -> dict[str, str]:
    attrs = {}
    for item in value.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
        elif " " in item:
            key, val = item.split(" ", 1)
            val = val.strip('"')
        else:
            continue
        attrs[key] = val
    return attrs


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def longest_stop_free_orf_aa(seq: str, strand: str) -> int:
    nt = revcomp(seq) if strand == "-" else seq
    protein = translate_nt(nt)
    best = current = 0
    for aa in protein:
        if aa == "*":
            best = max(best, current)
            current = 0
        else:
            current += 1
    return max(best, current)


def query_span_from_hsp(genome: dict[str, str], contig: str, strand: str, hsp: dict, qlen: int):
    seq = genome.get(contig, "")
    if strand == "+":
        start = hsp["sstart"] - (hsp["qstart"] - 1) * 3
        end = start + qlen * 3 - 1
    else:
        end = hsp["sstart"] + (hsp["qstart"] - 1) * 3
        start = end - qlen * 3 + 1
    if not seq or start < 1 or end > len(seq):
        return "", start, end
    nt = seq[start - 1:end]
    return (revcomp(nt) if strand == "-" else nt), start, end


def terminal_stop_after_span(genome: dict[str, str], contig: str, strand: str, start: int, end: int) -> str:
    seq = genome.get(contig, "")
    if not seq or not start or not end:
        return "na"
    if strand == "+":
        if end + 3 > len(seq):
            return "no"
        codon = seq[end:end + 3].upper()
    else:
        if start <= 3:
            return "no"
        codon = revcomp(seq[start - 4:start - 1]).upper()
    return "yes" if codon in STOP_CODONS else "no"


def parse_transcripts(gff_path: Path) -> dict[str, dict]:
    transcripts = {}
    cds_counts = Counter()
    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_attrs(cols[8])
            if cols[2] == "CDS":
                for parent in attrs.get("Parent", "").split(","):
                    if parent:
                        cds_counts[parent] += 1
                continue
            if cols[2] not in {"mRNA", "transcript"}:
                continue
            tid = attrs.get("ID") or attrs.get("transcript_id")
            if not tid:
                continue
            transcripts[tid] = {
                "contig": cols[0],
                "source": cols[1],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "gene": attrs.get("Parent") or attrs.get("gene_id") or tid,
            }
    for tid, tx in transcripts.items():
        tx["cds_count"] = cds_counts.get(tid, 0)
    return transcripts


def load_query_data(outdir: Path):
    query_seqs = read_fasta(outdir / "00_query" / "benchmark_effectors.fa")
    query_lengths = {}
    with (outdir / "00_query" / "benchmark_effectors.length.tsv").open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row:
                query_lengths[row[0]] = int(row[1])
    return query_seqs, query_lengths, list(query_lengths)


def classify_genomic_locus(config: Config, genome: dict[str, str], query: str, qlen: int, cluster: dict, hsps: list[dict]):
    contig = cluster["contig"]
    strand = cluster["strand"]
    locus_seq = genome.get(contig, "")[cluster["start"] - 1:cluster["end"]]
    genomic_orf_aa = longest_stop_free_orf_aa(locus_seq, strand) if locus_seq else 0
    max_qcov = max(h["qcov"] for h in hsps)
    best_identity = max(h["pident"] for h in hsps)
    best_span = {
        "same_frame_query_identity": 0.0,
        "query_span_start": "",
        "query_span_end": "",
        "query_span_internal_stop": "na",
        "query_span_starts_with_m": "na",
        "query_span_terminal_stop": "na",
    }
    for hsp in hsps:
        if hsp["pident"] < config.hard_effector_pident_min or hsp["qstart"] != 1:
            continue
        nt, span_start, span_end = query_span_from_hsp(genome, contig, strand, hsp, qlen)
        if not nt or len(nt) < qlen * 3:
            continue
        protein = translate_nt(nt)
        candidate = {
            "same_frame_query_identity": aa_identity(query, protein),
            "query_span_start": span_start,
            "query_span_end": span_end,
            "query_span_internal_stop": "yes" if "*" in protein else "no",
            "query_span_starts_with_m": "yes" if protein.startswith("M") else "no",
            "query_span_terminal_stop": terminal_stop_after_span(genome, contig, strand, span_start, span_end),
        }
        if candidate["same_frame_query_identity"] > best_span["same_frame_query_identity"]:
            best_span = candidate

    hard = (
        max_qcov >= config.hard_genome_qcov_min
        and best_identity >= config.hard_effector_pident_min
        and best_span["same_frame_query_identity"] >= config.hard_effector_pident_min
        and best_span["query_span_internal_stop"] == "no"
        and best_span["query_span_starts_with_m"] == "yes"
        and best_span["query_span_terminal_stop"] == "yes"
    )
    if genomic_orf_aa < config.genomic_orf_aa_min:
        status, reason = "genomic_short_orf", "short_orf"
    elif hard:
        status, reason = "genomic_complete", ""
    elif best_identity >= config.hard_effector_pident_min and max_qcov >= config.hard_genome_qcov_min and best_span["query_span_starts_with_m"] == "no":
        status, reason = "genomic_start_codon_missing", "start_codon_missing"
    elif best_identity >= config.hard_effector_pident_min and max_qcov >= config.hard_genome_qcov_min and best_span["query_span_terminal_stop"] == "no":
        status, reason = "genomic_terminal_stop_missing", "terminal_stop_missing"
    elif best_identity >= config.hard_effector_pident_min and max_qcov >= config.genome_complete_qcov_min:
        status, reason = "genomic_shifted_or_frameshift", "start_stop_shift_or_frameshift"
    elif max_qcov >= config.genome_complete_qcov_min:
        status, reason = "genomic_low_identity", "low_identity"
    else:
        status, reason = "genomic_partial", "partial_tblastn"
    best_span.update({"hard_intact_locus": "yes" if status == "genomic_complete" else "no", "exclusion_reason": reason})
    return status, genomic_orf_aa, best_span


def build_loci(config: Config, query_seqs: dict[str, str], query_lengths: dict[str, int]) -> list[dict]:
    loci = []
    # Each isolate is visited exactly once, so only the current assembly is held
    # in memory. Caching all assemblies would need tens of gigabytes for a
    # pan-genome-scale run.
    genome: dict[str, str] = {}
    for path in sorted((config.outdir / "10_tblastn_genome").glob("*/*.tblastn.tsv")):
        if path.name.startswith("._"):
            continue
        isolate = path.name.replace(".tblastn.tsv", "")
        genome = read_fasta(config.assembly_dir / f"{isolate}.fa")
        grouped = defaultdict(list)
        with path.open(newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 13 or row[0] not in query_lengths:
                    continue
                pident = float(row[2])
                qcov = float(row[12])
                if qcov < config.genome_qcov_min or pident < config.effector_pident_min:
                    continue
                sstart, send = int(row[8]), int(row[9])
                strand = "+" if sstart <= send else "-"
                grouped[(row[0], row[1], strand)].append({
                    "start": min(sstart, send),
                    "end": max(sstart, send),
                    "sstart": sstart,
                    "send": send,
                    "qstart": int(row[6]),
                    "qend": int(row[7]),
                    "pident": pident,
                    "qcov": qcov,
                    "bitscore": float(row[11]),
                })
        for (effector, contig, strand), hsps in grouped.items():
            hsps.sort(key=lambda item: (item["start"], item["end"]))
            clusters = []
            for hsp in hsps:
                if not clusters or hsp["start"] > clusters[-1]["end"] + config.locus_merge_distance:
                    clusters.append({"contig": contig, "strand": strand, "start": hsp["start"], "end": hsp["end"], "hsps": [hsp]})
                else:
                    clusters[-1]["end"] = max(clusters[-1]["end"], hsp["end"])
                    clusters[-1]["hsps"].append(hsp)
            for cluster in clusters:
                status, orf_aa, strict = classify_genomic_locus(
                    config, genome, query_seqs[effector], query_lengths[effector], cluster, cluster["hsps"]
                )
                loci.append({
                    "isolate": isolate,
                    "effector": effector,
                    "contig": contig,
                    "start": cluster["start"],
                    "end": cluster["end"],
                    "strand": strand,
                    "tblastn_identity": max(h["pident"] for h in cluster["hsps"]),
                    "tblastn_qcov": max(h["qcov"] for h in cluster["hsps"]),
                    "tblastn_bitscore": max(h["bitscore"] for h in cluster["hsps"]),
                    "tblastn_hsp_count": len(cluster["hsps"]),
                    "genomic_orf_aa": orf_aa,
                    "genomic_status": status,
                    **strict,
                })
    loci.sort(key=lambda x: (x["isolate"], x["effector"], x["contig"], x["start"], x["end"]))
    copy_counters = Counter()
    for locus in loci:
        key = (locus["isolate"], locus["effector"])
        copy_counters[key] += 1
        locus["genomic_copy_number"] = copy_counters[key]
        locus["locus_id"] = f'{locus["effector"]}_copy_{copy_counters[key]}'
    return loci


def load_blastp_hits(config: Config) -> dict[tuple[str, str], list[dict]]:
    """Load the per-isolate BLASTP hits against the complete merged proteome."""
    hits = defaultdict(list)
    for path in sorted((config.outdir / "30_blastp").glob("*/*.blastp.tsv")):
        if path.name.startswith("._"):
            continue
        isolate = path.name.replace(".blastp.tsv", "")
        with path.open(newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 16:
                    continue
                qlen, slen = int(row[14]), int(row[15])
                hits[(isolate, row[0])].append({
                    "isolate": isolate,
                    "effector": row[0],
                    "subject": row[1],
                    "tid": row[1],
                    "pident": float(row[2]),
                    "qcov": float(row[12]),
                    "scov": float(row[13]),
                    "qlen": qlen,
                    "slen": slen,
                    "length_ratio": slen / qlen if qlen else 0.0,
                    "bitscore": float(row[11]),
                })
    return hits


def status_for_hit(config: Config, hit: dict | None) -> str:
    if hit is None:
        return "missed"
    if hit["qcov"] >= config.model_qcov_min and hit["length_ratio"] >= config.model_len_ratio_min and hit["pident"] >= config.hard_effector_pident_min:
        return "complete"
    return "partial"


def best_stage_hit(config: Config, locus: dict, stage: str, blastp_hits: dict, transcripts: dict) -> dict | None:
    candidates = []
    for hit in blastp_hits.get((locus["isolate"], locus["effector"]), []):
        tx = transcripts.get(locus["isolate"], {}).get(hit["tid"])
        if not tx or tx["source"] not in STAGE_SOURCES[stage]:
            continue
        if tx["contig"] != locus["contig"] or tx["strand"] != locus["strand"]:
            continue
        ov = overlap_len(locus["start"], locus["end"], tx["start"], tx["end"])
        if ov <= 0:
            continue
        candidate = dict(hit)
        candidate.update({
            "gene": tx["gene"],
            "transcript": hit["tid"],
            "model_contig": tx["contig"],
            "model_source": tx["source"],
            "model_start": tx["start"],
            "model_end": tx["end"],
            "model_strand": tx["strand"],
            "model_cds_count": tx["cds_count"],
            "overextended_model": "yes" if hit["length_ratio"] > config.overextended_len_ratio else "no",
            "overlap_bp": ov,
        })
        candidates.append(candidate)
    if not candidates:
        return None
    return max(candidates, key=lambda h: (h["qcov"], -abs(1.0 - h["length_ratio"]), h["bitscore"], h["overlap_bp"]))


def hit_fields(hit: dict | None) -> list:
    if hit is None:
        return ["", "", "", "", "", "", "", 0, 0, 0, 0, "", 0, 0, ""]
    return [
        hit["gene"], hit["transcript"], hit["model_contig"], hit["model_source"],
        hit["model_start"], hit["model_end"], hit["model_strand"], f'{hit["pident"]:.2f}',
        f'{hit["qcov"]:.2f}', f'{hit["length_ratio"]:.3f}', f'{hit["bitscore"]:.1f}',
        hit["subject"], hit["overlap_bp"], hit["model_cds_count"], hit["overextended_model"],
    ]


def rescue_category(locus: dict, statuses: dict[str, str]) -> str:
    if locus["hard_intact_locus"] != "yes":
        return {
            "genomic_partial": "Genomic partial",
            "genomic_short_orf": "Genomic short ORF",
            "genomic_shifted_or_frameshift": "Genomic shifted/frameshift",
            "genomic_start_codon_missing": "Genomic start codon missing",
            "genomic_terminal_stop_missing": "Genomic terminal stop missing",
            "genomic_low_identity": "Genomic low identity",
        }.get(locus["genomic_status"], "Excluded genomic locus")
    if statuses["AUGUSTUS"] == "complete":
        return "AUGUSTUS"
    if statuses["AUGUSTUS_Helixer"] == "complete":
        return "Helixer rescue"
    if statuses["Final"] == "complete":
        return "miniprot rescue"
    return "Still missed"


def write_outputs(config: Config, loci: list[dict], effectors: list[str]) -> None:
    summary_dir = config.outdir / "40_summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    transcripts = {
        gff.name.replace("_merged.gff", ""): parse_transcripts(gff)
        for gff in sorted(config.merged_gff_dir.glob("*_merged.gff"))
        if not gff.name.startswith("._")
    }
    blastp_hits = load_blastp_hits(config)
    locus_path = summary_dir / "effector_locus_level_master.tsv"
    header = [
        "isolate", "effector", "locus_id", "contig", "start", "end", "strand",
        "genomic_copy_number", "genomic_status", "tblastn_identity", "tblastn_qcov",
        "tblastn_bitscore", "tblastn_hsp_count", "genomic_orf_aa",
        "same_frame_query_identity", "query_span_start", "query_span_end",
        "query_span_internal_stop", "query_span_starts_with_m", "query_span_terminal_stop",
        "hard_intact_locus", "exclusion_reason",
    ]
    for stage in STAGES:
        header.extend([
            f"{stage}_gene", f"{stage}_transcript", f"{stage}_model_contig",
            f"{stage}_model_source", f"{stage}_model_start", f"{stage}_model_end",
            f"{stage}_model_strand", f"{stage}_pident", f"{stage}_qcov",
            f"{stage}_length_ratio", f"{stage}_bitscore", f"{stage}_subject",
            f"{stage}_overlap_bp", f"{stage}_cds_count", f"{stage}_overextended_model",
            f"{stage}_status",
        ])
    header.append("rescue_method")

    rows = []
    for locus in loci:
        row = [
            locus["isolate"], locus["effector"], locus["locus_id"], locus["contig"],
            locus["start"], locus["end"], locus["strand"], locus["genomic_copy_number"],
            locus["genomic_status"], f'{locus["tblastn_identity"]:.2f}',
            f'{locus["tblastn_qcov"]:.2f}', f'{locus["tblastn_bitscore"]:.1f}',
            locus["tblastn_hsp_count"], locus["genomic_orf_aa"],
            f'{locus["same_frame_query_identity"]:.2f}', locus["query_span_start"],
            locus["query_span_end"], locus["query_span_internal_stop"],
            locus["query_span_starts_with_m"], locus["query_span_terminal_stop"],
            locus["hard_intact_locus"], locus["exclusion_reason"],
        ]
        statuses = {}
        for stage in STAGES:
            hit = best_stage_hit(config, locus, stage, blastp_hits, transcripts)
            statuses[stage] = status_for_hit(config, hit)
            row.extend(hit_fields(hit))
            row.append(statuses[stage])
        row.append(rescue_category(locus, statuses))
        rows.append(row)

    with locus_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Wrote {locus_path}")

    dict_rows = [dict(zip(header, row)) for row in rows]
    write_copy_summary(summary_dir, dict_rows, effectors)
    write_count_tables(summary_dir, dict_rows, effectors)


def write_copy_summary(summary_dir: Path, rows: list[dict], effectors: list[str]) -> None:
    isolates = sorted({row["isolate"] for row in rows})
    summary = defaultdict(Counter)
    for row in rows:
        key = (row["isolate"], row["effector"])
        summary[key]["genome_n"] += 1
        summary[key][f'{row["genomic_status"]}_n'] += 1
        if row["hard_intact_locus"] != "yes":
            summary[key]["excluded_genomic_locus_n"] += 1
            continue
        summary[key]["hard_intact_locus_n"] += 1
        for stage in STAGES:
            summary[key][f"{stage}_complete_n"] += row[f"{stage}_status"] == "complete"
        summary[key]["Helixer_rescue_n"] += row["rescue_method"] == "Helixer rescue"
        summary[key]["miniprot_rescue_n"] += row["rescue_method"] == "miniprot rescue"
        summary[key]["still_missed_n"] += row["rescue_method"] == "Still missed"

    fields = [
        "isolate", "effector", "genome_n", "hard_intact_locus_n", "excluded_genomic_locus_n",
        "genomic_complete_n", "genomic_partial_n", "genomic_short_orf_n",
        "genomic_shifted_or_frameshift_n", "genomic_start_codon_missing_n",
        "genomic_terminal_stop_missing_n", "genomic_low_identity_n",
        "AUGUSTUS_complete_n", "AUGUSTUS_Helixer_complete_n", "Final_complete_n",
        "Helixer_rescue_n", "miniprot_rescue_n", "still_missed_n",
    ]
    path = summary_dir / "effector_isolate_copy_summary.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        for isolate in isolates:
            for effector in effectors:
                values = summary.get((isolate, effector), {})
                writer.writerow([isolate, effector, *[values.get(field, 0) for field in fields[2:]]])
    print(f"Wrote {path}")


def write_count_tables(summary_dir: Path, rows: list[dict], effectors: list[str]) -> None:
    categories = [
        "AUGUSTUS", "Helixer rescue", "miniprot rescue", "Still missed",
        "Genomic partial", "Genomic short ORF", "Genomic shifted/frameshift",
        "Genomic start codon missing", "Genomic terminal stop missing",
        "Genomic low identity", "Excluded genomic locus",
    ]
    counts = defaultdict(Counter)
    stage_counts = defaultdict(Counter)
    for row in rows:
        counts[row["effector"]][row["rescue_method"]] += 1
        if row["hard_intact_locus"] == "yes":
            stage_counts[(row["effector"], "Genome")]["complete"] += 1
            for stage in STAGES:
                stage_counts[(row["effector"], stage)][row[f"{stage}_status"]] += 1
        else:
            stage_counts[(row["effector"], "Genome")]["excluded"] += 1

    path = summary_dir / "effector_locus_rescue_counts_by_effector.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["effector", "genomic_loci", *categories])
        for effector in effectors:
            total = sum(counts[effector].values())
            writer.writerow([effector, total, *[counts[effector].get(c, 0) for c in categories]])
    print(f"Wrote {path}")

    with (summary_dir / "effector_locus_stage_status_counts_by_effector.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["effector", "stage", "complete", "partial", "missed", "excluded"])
        for effector in effectors:
            for stage in ("Genome", *STAGES):
                values = stage_counts[(effector, stage)]
                writer.writerow([effector, stage, values["complete"], values["partial"], values["missed"], values["excluded"]])
    print(f"Wrote {summary_dir / 'effector_locus_stage_status_counts_by_effector.tsv'}")


def parse_args() -> Config:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    paths = parser.add_argument_group("input and output paths")
    paths.add_argument("--outdir", required=True, type=Path,
                       help="benchmark output directory written by benchmark.sh")
    paths.add_argument("--assembly-dir", required=True, type=Path,
                       help="directory of genome assemblies named <isolate>.fa")
    paths.add_argument("--merged-gff-dir", required=True, type=Path,
                       help="directory of merged annotations named <isolate>_merged.gff")

    screen = parser.add_argument_group("candidate genomic hit screening")
    screen.add_argument("--genome-qcov-min", type=float, default=50,
                        help="minimum tblastn query coverage (%%) to keep a candidate HSP")
    screen.add_argument("--effector-pident-min", type=float, default=85,
                        help="minimum tblastn amino acid identity (%%) to keep a candidate HSP")
    screen.add_argument("--locus-merge-distance", type=int, default=500,
                        help="maximum gap (bp) between HSPs merged into one genomic locus")

    intact = parser.add_argument_group("intact benchmark locus definition")
    intact.add_argument("--hard-genome-qcov-min", type=float, default=100,
                        help="query coverage (%%) required for an intact benchmark locus")
    intact.add_argument("--hard-effector-pident-min", type=float, default=95,
                        help="amino acid identity (%%) required for an intact benchmark locus")
    intact.add_argument("--genome-complete-qcov-min", type=float, default=90,
                        help="query coverage (%%) above which a non-intact locus is reported as shifted or low identity rather than partial")
    intact.add_argument("--genomic-orf-aa-min", type=int, default=50,
                        help="minimum stop-free ORF length (aa) at the locus")

    model = parser.add_argument_group("complete annotation definition")
    model.add_argument("--model-qcov-min", type=float, default=100,
                       help="minimum BLASTP query coverage (%%) of the gene model's protein")
    model.add_argument("--model-len-ratio-min", type=float, default=1.0,
                       help="minimum model/query protein length ratio")
    model.add_argument("--overextended-len-ratio", type=float, default=1.2,
                       help="model/query protein length ratio above which a model is flagged as overextended")

    return Config(**vars(parser.parse_args()))


def main() -> None:
    config = parse_args()
    query_seqs, query_lengths, effectors = load_query_data(config.outdir)
    loci = build_loci(config, query_seqs, query_lengths)
    write_outputs(config, loci, effectors)


if __name__ == "__main__":
    main()
