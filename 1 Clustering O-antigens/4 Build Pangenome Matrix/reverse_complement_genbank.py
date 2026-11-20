#!/usr/bin/env python3
"""Reverse complement every record in a multi-record GenBank file.

Reads a multi-record GenBank input file, reverse complements each record
(sequence, features and feature locations), and writes the results to a
new multi-record GenBank file. Record identifiers, descriptions,
annotations, database cross-references and qualifiers are preserved.

Usage:
    python reverse_complement_genbank.py input.gbff [output.gbff]

If no output path is given, writes <input_stem>_rc.gbff next to the input.
"""
import argparse
from pathlib import Path

from Bio import SeqIO


def reverse_complement_record(record):
    """Return a reverse-complemented copy of a SeqRecord, keeping metadata."""
    return record.reverse_complement(
        id=True,
        name=True,
        description=True,
        dbxrefs=True,
        annotations=True,
        letter_annotations=True,
        features=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Reverse complement all records in a multi-record GenBank file."
    )
    parser.add_argument("input", help="input multi-record GenBank file (.gbff/.gb/.gbk)")
    parser.add_argument(
        "output",
        nargs="?",
        help="output GenBank file (default: <input_stem>_rc.gbff beside the input)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        parser.error(f"input file not found: {input_path}")

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(f"{input_path.stem}_rc.gbff")

    records = [reverse_complement_record(record) for record in SeqIO.parse(input_path, "genbank")]
    if not records:
        parser.error(f"no GenBank records found in: {input_path}")

    count = SeqIO.write(records, output_path, "genbank")
    print(f"Reverse complemented {count} records: {input_path} -> {output_path}")


if __name__ == "__main__":
    main()
