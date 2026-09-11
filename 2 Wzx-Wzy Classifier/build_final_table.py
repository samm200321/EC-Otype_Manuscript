#!/usr/bin/env python3
"""
Build the final per-protein assignment table:
  - proteins matched to wzx/wzy at >=identity threshold -> that label
  - unmatched proteins -> "other_cluster_<N>" based on mmseqs easy-cluster output
    (proteins with no cluster file, e.g. because clustering was skipped, get "other")
"""
import argparse
import os
import pandas as pd
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unlabeled-fasta", required=True)
    ap.add_argument("--matched", required=True, help="matched_labels.tsv from assign_labels.py")
    ap.add_argument("--cluster-tsv", required=True,
                     help="mmseqs easy-cluster *_cluster.tsv output (rep<TAB>member)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    all_ids = [rec.id for rec in SeqIO.parse(args.unlabeled_fasta, "fasta")]
    df = pd.DataFrame({"id": all_ids})

    matched = pd.read_csv(args.matched, sep="\t")
    matched = matched.rename(columns={"query": "id"})[["id", "label", "pident"]]

    df = df.merge(matched, on="id", how="left")

    # Assign cluster-based labels to whatever wasn't matched to wzx/wzy.
    if os.path.exists(args.cluster_tsv) and os.path.getsize(args.cluster_tsv) > 0:
        clusters = pd.read_csv(args.cluster_tsv, sep="\t", names=["representative", "id"])
        # Give each distinct representative a stable, readable cluster name.
        rep_to_clusterid = {
            rep: f"other_cluster_{i+1}"
            for i, rep in enumerate(sorted(clusters["representative"].unique()))
        }
        clusters["cluster_label"] = clusters["representative"].map(rep_to_clusterid)
        df = df.merge(clusters[["id", "representative", "cluster_label"]], on="id", how="left")
    else:
        df["representative"] = None
        df["cluster_label"] = None

    # Final label: wzx/wzy where matched, else the cluster label, else "other_unclustered"
    # (this last case shouldn't normally occur unless clustering was skipped entirely).
    df["final_label"] = df["label"].fillna(df["cluster_label"]).fillna("other_unclustered")

    df = df[["id", "final_label", "label", "pident", "cluster_label", "representative"]]
    df.to_csv(args.out, sep="\t", index=False)

    print(df["final_label"].apply(
        lambda x: "wzx" if x == "wzx" else "wzy" if x == "wzy" else "other"
    ).value_counts())
    print(f"\nWrote {len(df)} rows to {args.out}")


if __name__ == "__main__":
    main()
