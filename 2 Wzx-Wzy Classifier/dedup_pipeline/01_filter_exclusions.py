#!/usr/bin/env python3
"""Step 1+2: filter data/all_prot_2.faa by exclusion lists.

Excludes:
  1) proteins listed in the "query" column of wzx_wzy_labels/matched_labels.tsv
  2) proteins listed in data/short_proteins_lt100aa.csv

Writes: dedup_pipeline/results/filtered_all_prot.faa
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ALL_PROT = ROOT / "data" / "all_prot_2.faa"
MATCHED_LABELS = ROOT / "wzx_wzy_labels" / "matched_labels.tsv"
SHORT_PROTS = ROOT / "data" / "short_proteins_lt100aa.csv"
OUT_DIR = ROOT / "dedup_pipeline" / "results"
OUT_FASTA = OUT_DIR / "filtered_all_prot.faa"


def load_query_exclusions(path: Path) -> set:
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["query"] for row in reader}


def load_short_exclusions(path: Path) -> set:
    with open(path) as fh:
        reader = csv.DictReader(fh)
        return {row["protein_name"] for row in reader}


def read_fasta(path: Path):
    """Yield (header_without_>, sequence) tuples."""
    header, seq_lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def main():
    query_excl = load_query_exclusions(MATCHED_LABELS)
    short_excl = load_short_exclusions(SHORT_PROTS)
    combined_excl = query_excl | short_excl

    n_total = n_kept = n_excl_query = n_excl_short = 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_FASTA, "w") as out:
        for header, seq in read_fasta(ALL_PROT):
            n_total += 1
            name = header.split()[0]
            if name in combined_excl:
                if name in query_excl:
                    n_excl_query += 1
                if name in short_excl:
                    n_excl_short += 1
                continue
            n_kept += 1
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

    print(f"Input proteins:            {n_total}")
    print(f"Excluded (matched query):  {n_excl_query}")
    print(f"Excluded (short <100 aa):  {n_excl_short}")
    print(f"Excluded (union):          {n_total - n_kept}")
    print(f"Kept:                      {n_kept}")
    print(f"Wrote {OUT_FASTA}")


if __name__ == "__main__":
    main()
