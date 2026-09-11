#!/usr/bin/env python3
"""
Parse mmseqs search hits (query vs. labeled wzx/wzy reference set), keep the
best-identity hit per query, and assign wzx/wzy labels to queries that clear
the identity + coverage thresholds. Queries with no qualifying hit are
written out to a separate FASTA for downstream clustering.

Assumes labeled FASTA headers encode the label, e.g.:
    >wzx_Ecoli_O157
    >wzy_Klebsiella_K1
    >wzx|some_other_id

Adjust `infer_label_from_id()` below if your labeled FASTA is structured
differently (e.g. a separate id->label mapping TSV).
"""
import argparse
import re
import pandas as pd
from Bio import SeqIO


def infer_label_from_id(seq_id: str) -> str:
    """Infer 'wzx' or 'wzy' from a labeled-set sequence ID. Customize as needed."""
    lower = seq_id.lower()
    if re.search(r"\bwzx\b", lower) or lower.startswith("wzx"):
        return "wzx"
    if re.search(r"\bwzy\b", lower) or lower.startswith("wzy"):
        return "wzy"
    raise ValueError(
        f"Could not infer wzx/wzy label from labeled sequence ID: {seq_id!r}. "
        "Edit infer_label_from_id() to match your header format, or supply a "
        "separate id->label mapping."
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="mmseqs convertalis output TSV")
    ap.add_argument("--labeled-fasta", required=True)
    ap.add_argument("--unlabeled-fasta", required=True)
    ap.add_argument("--id-threshold", type=float, default=0.95,
                     help="Minimum fractional identity (0-1) to accept a match")
    ap.add_argument("--cov-threshold", type=float, default=0.8,
                     help="Minimum coverage (0-1) to accept a match")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    cols = ["query", "target", "pident", "alnlen", "qcov", "tcov", "evalue", "bits"]
    hits = pd.read_csv(args.hits, sep="\t", names=cols)

    # mmseqs convertalis reports pident as a fraction (0-1) by default with
    # --format-output; if your build reports 0-100 instead, divide by 100 here.
    # Sanity check:
    if hits["pident"].max() > 1.5:
        print("Note: pident appears to be on a 0-100 scale; rescaling to 0-1.")
        hits["pident"] = hits["pident"] / 100.0

    # Apply thresholds, then keep the single best-identity hit per query.
    passing = hits[
        (hits["pident"] >= args.id_threshold) &
        (hits[["qcov", "tcov"]].min(axis=1) >= args.cov_threshold)
    ].copy()

    best_hits = (
        passing.sort_values("pident", ascending=False)
        .drop_duplicates(subset="query", keep="first")
    )

    # Map matched labeled-set target IDs -> wzx/wzy label.
    labeled_ids = [rec.id for rec in SeqIO.parse(args.labeled_fasta, "fasta")]
    label_map = {seq_id: infer_label_from_id(seq_id) for seq_id in labeled_ids}

    best_hits["label"] = best_hits["target"].map(label_map)
    unmapped_targets = best_hits[best_hits["label"].isna()]
    if len(unmapped_targets):
        raise RuntimeError(
            f"{len(unmapped_targets)} hits matched a labeled-set target ID that "
            "could not be mapped to wzx/wzy. Check infer_label_from_id()."
        )

    matched_out = best_hits[["query", "target", "pident", "qcov", "tcov", "label"]]
    matched_out.to_csv(f"{args.outdir}/matched_labels.tsv", sep="\t", index=False)
    print(f"Matched {len(matched_out)} queries to wzx/wzy "
          f"(wzx={sum(matched_out.label == 'wzx')}, wzy={sum(matched_out.label == 'wzy')})")

    # Everything else -> unmatched, write out for clustering.
    all_unlabeled_ids = {rec.id for rec in SeqIO.parse(args.unlabeled_fasta, "fasta")}
    matched_ids = set(matched_out["query"])
    unmatched_ids = all_unlabeled_ids - matched_ids
    print(f"{len(unmatched_ids)} / {len(all_unlabeled_ids)} queries unmatched -> clustering separately")

    with open(f"{args.outdir}/unmatched.fasta", "w") as out_f:
        for rec in SeqIO.parse(args.unlabeled_fasta, "fasta"):
            if rec.id in unmatched_ids:
                SeqIO.write(rec, out_f, "fasta")


if __name__ == "__main__":
    main()
