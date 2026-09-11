#!/usr/bin/env python3
"""Build an exclusion list of locus tags from a GenBank file.

A CDS is excluded if:
  - its /product contains "polymerase" or "flippase" (case-insensitive), or
  - it has no /product qualifier.
"""

import csv

from Bio import SeqIO

GBFF_FILE = "data/EC_Oag_v0.2.gbff"
OUTPUT_FILE = "exclusion_list.csv"
KEYWORDS = ("polymerase", "flippase", "Wzy", "Wzx", "epsg")


def main():
    excluded = []
    seen = set()
    total_cds = 0

    for record in SeqIO.parse(GBFF_FILE, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            total_cds += 1

            locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
            products = feature.qualifiers.get("product")

            product = products[0].lower() if products else None
            exclude = product is None or any(kw in product for kw in KEYWORDS)

            if product == "wzy":
                exclude = True
            elif product == "wzx":
                exclude = True

            if exclude and locus_tag and locus_tag not in seen:
                seen.add(locus_tag)
                excluded.append(locus_tag)
                print(product)

    with open(OUTPUT_FILE, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["locus_tag"])
        writer.writerows([[tag] for tag in excluded])

    print(f"Scanned {total_cds} CDS features across {len(list(SeqIO.parse(GBFF_FILE, 'genbank')))} records.")
    print(f"Wrote {len(excluded)} locus tags to {OUTPUT_FILE}.")


if __name__ == "__main__":
    main()
