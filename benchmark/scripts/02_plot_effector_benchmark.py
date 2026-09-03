#!/usr/bin/env python3
"""Plot copy-aware effector annotation benchmark results.

Reads the summary tables written by 01_build_locus_tables.py and writes the
benchmark figures.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# The keys are the column prefixes written by 01_build_locus_tables.py, which
# follow the GFF source values; the titles use BRAKER, the program that writes
# `AUGUSTUS` into that column.
STAGES = [
    ("Genomic copies", "genome"),
    ("Copies annotated by BRAKER", "AUGUSTUS"),
    ("Copies annotated by BRAKER + Helixer", "AUGUSTUS_Helixer"),
    ("Copies annotated by BRAKER + Helixer + miniprot", "Final"),
]

# The rescue_method values written by 01_build_locus_tables.py are kept stable in
# the tables; only their display labels are set here.
RESCUE_LABELS = {
    "AUGUSTUS": "Annotated by BRAKER",
    "Helixer rescue": "Added by Helixer",
    "miniprot rescue": "Added by miniprot",
    "Still missed": "Not annotated",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_all(fig, plot_dir: Path, stem: str, formats: list[str], dpi: int) -> list[Path]:
    paths = []
    for ext in formats:
        path = plot_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def plot_copy_aware_counts(rows: list[dict[str, str]], plot_dir: Path, formats: list[str], dpi: int) -> list[Path]:
    effectors = list(dict.fromkeys(row["effector"] for row in rows))
    isolates = sorted({row["isolate"] for row in rows})
    complete_by_key = {}
    excluded_by_key = defaultdict(int)
    for effector in effectors:
        for isolate in isolates:
            loci = [row for row in rows if row["effector"] == effector and row["isolate"] == isolate]
            intact = [row for row in loci if row["hard_intact_locus"] == "yes"]
            complete_by_key[(effector, isolate, "genome")] = len(intact)
            excluded_by_key[(effector, isolate)] = sum(row["hard_intact_locus"] != "yes" for row in loci)
            for stage in ("AUGUSTUS", "AUGUSTUS_Helixer", "Final"):
                complete_by_key[(effector, isolate, stage)] = sum(row[f"{stage}_status"] == "complete" for row in intact)

    fig_w = max(9, 0.17 * len(isolates) + 2.6)
    fig_h = max(8, 1.55 * len(STAGES) + 0.28 * len(effectors))
    fig, axes = plt.subplots(len(STAGES), 1, figsize=(fig_w, fig_h), sharex=True)
    present_color = "#2ca58d"
    incomplete_color = "#e34a33"
    absent_color = "#d9d9d9"

    for ax, (title, stage) in zip(axes, STAGES):
        ax.set_xlim(-0.5, len(isolates) - 0.5)
        ax.set_ylim(len(effectors) - 0.5, -0.5)
        ax.set_title(title, loc="left", pad=6)
        ax.set_ylabel("Benchmark effector")
        ax.set_yticks(range(len(effectors)))
        ax.set_yticklabels(effectors, fontsize=8)
        ax.set_xticks([x - 0.5 for x in range(1, len(isolates))], minor=True)
        ax.set_yticks([y - 0.5 for y in range(1, len(effectors))], minor=True)
        ax.grid(which="minor", color="#d9d9d9", linewidth=0.35)
        ax.tick_params(which="minor", bottom=False, left=False)
        for y, effector in enumerate(effectors):
            for x, isolate in enumerate(isolates):
                genome_n = complete_by_key.get((effector, isolate, "genome"), 0)
                value = complete_by_key.get((effector, isolate, stage), 0)
                if stage == "genome":
                    color = present_color if genome_n > 0 else absent_color
                    label = str(genome_n) if genome_n > 0 else ""
                else:
                    if genome_n == 0:
                        color = absent_color
                        label = ""
                    elif value == genome_n:
                        color = present_color
                        label = str(value)
                    else:
                        color = incomplete_color
                        label = str(value)
                ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=color, edgecolor="none"))
                if label:
                    ax.text(x, y, label, ha="center", va="center", fontsize=5.5, color="#111111")

    axes[-1].set_xticks(range(len(isolates)))
    axes[-1].set_xticklabels(isolates, rotation=60, ha="right", fontsize=5.5)
    axes[-1].set_xlabel("Isolate", labelpad=8)
    legend_handles = [
        Patch(facecolor=present_color, edgecolor="none", label="all copies annotated"),
        Patch(facecolor=incomplete_color, edgecolor="none", label="one or more copies unannotated"),
        Patch(facecolor=absent_color, edgecolor="none", label="no intact copy detected"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, 0.015), ncol=3, frameon=False)
    fig.subplots_adjust(left=0.08, right=0.92, top=0.97, bottom=0.20, hspace=0.42)
    return save_all(fig, plot_dir, "effector_annotation_count_heatmap", formats, dpi)


def plot_rescue_counts(rows: list[dict[str, str]], plot_dir: Path, formats: list[str], dpi: int) -> list[Path]:
    effectors = list(dict.fromkeys(row["effector"] for row in rows))
    categories = ["AUGUSTUS", "Helixer rescue", "miniprot rescue", "Still missed"]
    colors = {
        "AUGUSTUS": "#2ca58d",
        "Helixer rescue": "#4c78a8",
        "miniprot rescue": "#8e6bbe",
        "Still missed": "#e34a33",
    }
    counts = defaultdict(Counter)
    for row in rows:
        if row["hard_intact_locus"] == "yes":
            counts[row["effector"]][row["rescue_method"]] += 1

    fig_h = max(4, 0.35 * len(effectors) + 1.8)
    fig, ax = plt.subplots(figsize=(8, fig_h))
    left = [0] * len(effectors)
    y_positions = range(len(effectors))
    for category in categories:
        widths = [counts[effector][category] for effector in effectors]
        ax.barh(y_positions, widths, left=left, color=colors[category], edgecolor="white", linewidth=0.8,
                label=RESCUE_LABELS[category])
        for y, lft, width in zip(y_positions, left, widths):
            if width >= 3:
                ax.text(lft + width / 2, y, str(width), ha="center", va="center", fontsize=7)
        left = [lft + width for lft, width in zip(left, widths)]
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(effectors)
    ax.invert_yaxis()
    ax.set_xlabel("Genomic copies")
    ax.set_ylabel("Benchmark effector")
    ax.set_title("Source set at which each genomic copy first becomes annotated")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=4)
    fig.tight_layout()
    return save_all(fig, plot_dir, "effector_rescue_counts_by_effector", formats, dpi)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--summary-dir", required=True, type=Path,
                        help="directory holding effector_locus_level_master.tsv")
    parser.add_argument("--plot-dir", required=True, type=Path,
                        help="directory to write the figures into")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"],
                        help="figure file formats to write")
    parser.add_argument("--dpi", type=int, default=300,
                        help="resolution of the raster output")
    args = parser.parse_args()
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    table = args.summary_dir / "effector_locus_level_master.tsv"
    if not table.exists():
        raise SystemExit(f"ERROR: summary table not found: {table}")
    rows = read_tsv(table)
    paths = []
    paths.extend(plot_copy_aware_counts(rows, args.plot_dir, args.formats, args.dpi))
    paths.extend(plot_rescue_counts(rows, args.plot_dir, args.formats, args.dpi))
    for path in paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
