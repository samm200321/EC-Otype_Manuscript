#!/usr/bin/env python3
"""Step 3b: histogram of pairwise identities.

Reads dedup_pipeline/results/pairwise_identities.tsv and plots the
distribution of pairwise % identities. Candidate de-duplication
thresholds are marked so one can be chosen visually; the selected
threshold is passed to 03_cluster_dedup.py.

Writes: dedup_pipeline/results/pairwise_identity_histogram.png
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "dedup_pipeline" / "results"
IN_TSV = RESULTS / "pairwise_identities.tsv"
OUT_PNG = RESULTS / "pairwise_identity_histogram.png"

CANDIDATE_THRESHOLDS = [0.5, 0.7, 0.8, 0.9, 0.95]


def main():
    pidents = []
    with open(IN_TSV) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pidents.append(float(row["pident"]))
    print(f"Loaded {len(pidents)} pairwise identities from {IN_TSV}")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(pidents, bins=100, color="steelblue", edgecolor="none")
    ax.set_xlabel("Pairwise sequence identity")
    ax.set_ylabel("Count")
    ax.set_title(f"All-vs-all pairwise identities (n hits = {len(pidents)})")
    ax.set_yscale("log")

    for thr in CANDIDATE_THRESHOLDS:
        ax.axvline(thr, linestyle="--", linewidth=1,
                   label=f"{thr:.2f}")
    ax.legend(title="Candidate thresholds")

    # quantile summary
    pidents.sort()
    n = len(pidents)
    qs = {q: pidents[min(n - 1, int(q * n))] for q in (0.5, 0.9, 0.99)}
    txt = "\n".join(f"p{int(q*100)}: {v:.3f}" for q, v in qs.items())
    ax.text(0.02, 0.95, txt, transform=ax.transAxes, va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    plt.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")
    print("Median/p90/p99 identity:",
          {f"p{int(q*100)}": round(v, 3) for q, v in qs.items()})


if __name__ == "__main__":
    main()
