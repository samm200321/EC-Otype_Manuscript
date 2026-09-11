#!/usr/bin/env python3
"""Filter a multi-record GenBank file to keep only records with >= N CDS features."""
import sys
from Bio import SeqIO


def main():
    input_gbff = "ISfinder_results/O-antigen_db.reverse_complement.is_removed.no_transposase.gbff"
    output_gbff = "ISfinder_results/O-antigen_db.reverse_complement.is_removed.no_transposase.min3orfs.gbff"
    min_orfs = 3

    kept = 0
    skipped = 0
    with open(output_gbff, "w") as out_handle:
        for record in SeqIO.parse(input_gbff, "genbank"):
            cds_count = sum(1 for f in record.features if f.type == "CDS")
            if cds_count >= min_orfs:
                SeqIO.write(record, out_handle, "genbank")
                kept += 1
            else:
                skipped += 1

    print(f"Kept {kept} records with >= {min_orfs} ORFs")
    print(f"Skipped {skipped} records with < {min_orfs} ORFs")
    print(f"Output: {output_gbff}")


if __name__ == "__main__":
    main()
