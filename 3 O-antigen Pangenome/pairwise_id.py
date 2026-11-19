#!/usr/bin/env python3
"""
pairwise_identity_mmseqs.py

Compute all-vs-all pairwise amino acid identities for the sequences in a
multi-FASTA file using MMseqs2, then plot the distribution as a histogram.

Requirements:
    - mmseqs2 installed and on PATH (https://github.com/soedinglab/MMseqs2)
    - Python packages: pandas, matplotlib, seaborn

Usage:
    python pairwise_identity_mmseqs.py input.fasta -o results \
        --sensitivity 7.5 --threads 4 --min-cov 0.5 --bins 50

Output:
    results/pairwise_hits.tsv   -> raw mmseqs all-vs-all hits table
    results/pairwise_identity_histogram.png -> histogram + KDE plot
    results/pairwise_identity_histogram.pdf -> same figure as a vector PDF with
        Calibri text embedded as editable Type 42 fonts (for Adobe Illustrator)
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
import seaborn as sns


MMSEQS_SEARCH_COLUMNS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]


def check_mmseqs():
    if shutil.which("mmseqs") is None:
        sys.exit(
            "ERROR: 'mmseqs' was not found on your PATH.\n"
            "Install MMseqs2, e.g.:\n"
            "  conda install -c bioconda mmseqs2\n"
            "or see https://github.com/soedinglab/MMseqs2 for other install options."
        )


def run(cmd, **kwargs):
    print(f"$ {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True, **kwargs)


def run_mmseqs_all_vs_all(fasta_path: Path, workdir: Path, threads: int,
                           sensitivity: float, min_cov: float) -> Path:
    """Run mmseqs createdb + search (all-vs-all) + convertalis. Returns path to TSV."""
    db = workdir / "seqDB"
    result_db = workdir / "resultDB"
    tmp_dir = workdir / "tmp"
    out_tsv = workdir / "pairwise_hits.tsv"
    tmp_dir.mkdir(exist_ok=True)

    run(["mmseqs", "createdb", str(fasta_path), str(db)])

    run([
        "mmseqs", "search", str(db), str(db), str(result_db), str(tmp_dir),
        "--threads", str(threads),
        "-s", str(sensitivity),
        "-c", str(min_cov),
        "--cov-mode", "0",
        # allow self hits so single-sequence edge cases don't error, filtered later
    ])

    run([
        "mmseqs", "convertalis", str(db), str(db), str(result_db), str(out_tsv),
        "--format-output", ",".join(MMSEQS_SEARCH_COLUMNS),
    ])

    return out_tsv


def load_and_filter(tsv_path: Path, drop_self: bool = True,
                     dedupe_pairs: bool = True) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", names=MMSEQS_SEARCH_COLUMNS)

    if drop_self:
        df = df[df["query"] != df["target"]]

    if dedupe_pairs:
        # MMseqs all-vs-all search reports both A->B and B->A (and these can
        # differ slightly since alignment is not symmetric). Keep one row per
        # unordered pair, taking the max identity seen for that pair.
        pair_key = df.apply(lambda r: tuple(sorted((r["query"], r["target"]))), axis=1)
        df = df.assign(_pair=pair_key)
        df = df.sort_values("pident", ascending=False).drop_duplicates("_pair")
        df = df.drop(columns="_pair")

    return df.reset_index(drop=True)


def configure_pdf_fonts():
    """Set matplotlib rcParams so the exported PDF opens in Adobe Illustrator
    with all text editable (Calibri, Type 42 TrueType embedding)."""
    # Type 42 (TrueType) embedding keeps text as real, editable text objects in
    # Illustrator; the default Type 3 embedding turns glyphs into un-editable shapes.
    plt.rcParams["pdf.fonttype"] = 42
    available_fonts = {f.name for f in font_manager.fontManager.ttflist}
    if "Calibri" in available_fonts:
        plt.rcParams["font.family"] = "Calibri"
    else:
        print("WARNING: 'Calibri' font not found on this system; "
              "falling back to the default sans-serif font in the PDF.")


def plot_histogram(df: pd.DataFrame, out_png: Path, bins: int):
    sns.set_theme(style="white")
    # Must be applied after set_theme, which resets font rcParams
    configure_pdf_fonts()
    plt.figure(figsize=(8.27, 5))
    ax = sns.histplot(
        data=df, x="pident", bins=bins, kde=True,
        color="#4C72B0", edgecolor="black", alpha=0.7,
        line_kws={"linewidth": 2, "color": "#C44E52"},
    )
    ax.set_xlabel("Pairwise amino acid identity (%)")
    ax.set_ylabel("Number of sequence pairs")
    #ax.set_title(f"Pairwise identity distribution (n={len(df)} pairs)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    out_pdf = out_png.with_suffix(".pdf")
    plt.savefig(out_pdf)
    print(f"Saved histogram to {out_png}")
    print(f"Saved Illustrator-editable PDF to {out_pdf}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", type=Path, help="Input multi-protein FASTA file")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("mmseqs_identity_out"),
                         help="Output directory (default: mmseqs_identity_out)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for mmseqs (default: 4)")
    parser.add_argument("--sensitivity", type=float, default=7.5,
                         help="mmseqs search sensitivity -s (default: 7.5, max sensitivity)")
    parser.add_argument("--min-cov", type=float, default=0.0,
                         help="Minimum alignment coverage 0-1 passed to mmseqs -c (default: 0, i.e. no filter)")
    parser.add_argument("--bins", type=int, default=50, help="Number of histogram bins (default: 50)")
    parser.add_argument("--keep-self-hits", action="store_true",
                         help="Keep self-vs-self hits (excluded by default)")
    parser.add_argument("--keep-tmp", action="store_true",
                         help="Keep MMseqs2 intermediate database files")
    args = parser.parse_args()

    if not args.fasta.exists():
        sys.exit(f"ERROR: input FASTA not found: {args.fasta}")

    check_mmseqs()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Use a temp working dir for mmseqs DB files, but keep final TSV/plot in outdir
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        raw_tsv = run_mmseqs_all_vs_all(
            fasta_path=args.fasta,
            workdir=workdir,
            threads=args.threads,
            sensitivity=args.sensitivity,
            min_cov=args.min_cov,
        )

        final_tsv = args.outdir / "pairwise_hits.tsv"
        shutil.copy(raw_tsv, final_tsv)

        if args.keep_tmp:
            kept_dir = args.outdir / "mmseqs_tmp"
            shutil.copytree(workdir, kept_dir, dirs_exist_ok=True)
            print(f"Kept mmseqs intermediate files in {kept_dir}")

    df = load_and_filter(final_tsv, drop_self=not args.keep_self_hits)

    if df.empty:
        sys.exit("No pairwise hits found (after filtering). "
                  "Try lowering --min-cov or check your input FASTA has >1 sequence.")

    stats_path = args.outdir / "pairwise_identity_stats.txt"
    with open(stats_path, "w") as fh:
        summary = df["pident"].describe()
        fh.write(str(summary))
    print(summary)

    hist_path = args.outdir / "pairwise_identity_histogram.png"
    plot_histogram(df, hist_path, bins=args.bins)

    filtered_tsv = args.outdir / "pairwise_identity_filtered.tsv"
    df.to_csv(filtered_tsv, sep="\t", index=False)
    print(f"Filtered/deduped pairwise table written to {filtered_tsv}")


if __name__ == "__main__":
    main()