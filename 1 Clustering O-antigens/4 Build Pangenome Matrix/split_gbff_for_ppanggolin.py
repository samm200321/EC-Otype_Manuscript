#!/usr/bin/env python3
"""Split a multi-record GenBank file into individual files and create a
PPanGGOLiN-compatible TSV manifest (accession<tab>filepath).
"""
from pathlib import Path
from Bio import SeqIO


def main():
    input_gbff = "min3/O-antigen_db.transposase_clean.min3orfs_rc.gbff"
    output_dir = Path("ppanggolin_input")
    output_dir.mkdir(exist_ok=True)

    records = []
    for record in SeqIO.parse(input_gbff, "genbank"):
        acc = record.id.replace("/", "_").replace(" ", "_")
        if not acc:
            acc = record.name.replace("/", "_").replace(" ", "_")
        out_path = output_dir / f"{acc}.gbff"
        SeqIO.write(record, out_path, "genbank")
        records.append((acc, str(out_path)))

    tsv_path = "ppanggolin_input.tsv"
    with open(tsv_path, "w") as fh:
        for acc, path in records:
            fh.write(f"{acc}\t{path}\n")

    print(f"Split {len(records)} records into {output_dir}/")
    print(f"Wrote manifest: {tsv_path}")
    print("First 5 lines:")
    with open(tsv_path) as fh:
        for line in fh.readlines()[:5]:
            print(" ", line.strip())


if __name__ == "__main__":
    main()
