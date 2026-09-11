#!/usr/bin/env python3
"""Step 3a: all-vs-all MMseqs2 search to obtain pairwise identities.

Combines the filtered all_prot set with the wzx/wzy reference proteins,
then runs a self-search to produce a table of pairwise % identities
(self-hits removed). The table feeds the histogram step and informs the
choice of clustering identity threshold.

Usage: 03_pairwise_identities.py [--min-seq-id 0.3]

Writes: dedup_pipeline/results/combined_input.faa
        dedup_pipeline/results/pairwise_identities.tsv
"""
import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "dedup_pipeline" / "results"
FILTERED = RESULTS / "filtered_all_prot.faa"
WZX = ROOT / "data" / "wzx_proteins.faa"
WZY = ROOT / "data" / "wzy_proteins.faa"
COMBINED = RESULTS / "combined_input.faa"
ALN_TSV = RESULTS / "pairwise_identities_raw.tsv"
OUT_TSV = RESULTS / "pairwise_identities.tsv"
TMP = RESULTS / "tmp_search"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-seq-id", type=float, default=0.3,
                    help="minimum identity to report in the search (default 0.3)")
    args = ap.parse_args()

    # concatenate filtered + wzx + wzy
    with open(COMBINED, "w") as out:
        for src in (FILTERED, WZX, WZY):
            out.write(src.read_text())
    n = sum(1 for line in open(COMBINED) if line.startswith(">"))
    print(f"Combined input: {n} sequences -> {COMBINED}")

    TMP.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs", "easy-search",
        str(COMBINED), str(COMBINED), str(ALN_TSV), str(TMP),
        "--min-seq-id", str(args.min_seq_id),
        "-s", "7.5",
        "--alignment-mode", "3",
        "--format-output", "query,target,pident,alnlen,qcov,tcov",
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # strip self-hits; scale pident (0-100) -> fraction
    n_hits = 0
    with open(ALN_TSV) as fh, open(OUT_TSV, "w") as out:
        out.write("query\ttarget\tpident\talnlen\tqcov\ttcov\n")
        for line in fh:
            q, t, pident, alnlen, qcov, tcov = line.rstrip("\n").split("\t")
            if q == t:
                continue
            pident = float(pident) / 100.0
            out.write(f"{q}\t{t}\t{pident:.4f}\t{alnlen}\t{qcov}\t{tcov}\n")
            n_hits += 1
    print(f"Pairwise hits (self-hits removed): {n_hits} -> {OUT_TSV}")


if __name__ == "__main__":
    main()
