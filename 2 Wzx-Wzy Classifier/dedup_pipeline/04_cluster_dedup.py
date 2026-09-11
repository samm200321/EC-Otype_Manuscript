#!/usr/bin/env python3
"""Step 3c: MMseqs2 clustering (de-duplication) at a chosen identity threshold.

Clusters the combined set (filtered all_prot + wzx + wzy) and writes the
representative sequences as the final de-duplicated multi-FASTA.

Usage: 04_cluster_dedup.py [--min-seq-id 0.9] [--cov 0.8]

Writes: dedup_pipeline/results/deduplicated_combined.faa
"""
import argparse
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "dedup_pipeline" / "results"
COMBINED = RESULTS / "combined_input.faa"
OUT_FASTA = RESULTS / "deduplicated_combined.faa"
PREFIX = RESULTS / "cluster"
TMP = RESULTS / "tmp_cluster"


def count_fasta(path: Path) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-seq-id", type=float, default=0.99,
                    help="clustering identity threshold (default 0.95)")
    ap.add_argument("--cov", type=float, default=1.00,
                    help="minimum alignment coverage (default 0.95)")
    args = ap.parse_args()

    TMP.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs", "easy-cluster",
        str(COMBINED), str(PREFIX), str(TMP),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.cov),
        "--cov-mode", "1",  # coverage of the shorter sequence
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    rep = Path(f"{PREFIX}_rep_seq.fasta")
    shutil.copy(rep, OUT_FASTA)

    n_in = count_fasta(COMBINED)
    n_out = count_fasta(OUT_FASTA)
    print(f"Input sequences:              {n_in}")
    print(f"Clusters (representatives):   {n_out}")
    print(f"Removed by de-duplication:    {n_in - n_out} "
          f"({100 * (n_in - n_out) / n_in:.2f}%)")
    print(f"Identity threshold:           {args.min_seq_id}")
    print(f"Wrote {OUT_FASTA}")


if __name__ == "__main__":
    main()
