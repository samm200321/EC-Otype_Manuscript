#!/usr/bin/env python3
"""Parse a GenBank file (e.g. Bakta .gbff) and physically delete the nucleotide
sequence of CDS features whose /product looks like a transposase / mobile
element protein, a phage / prophage protein, or an integrase / recombinase /
resolvase. Housekeeping recombination genes (RecA, XerC/XerD, RecFOR) are
protected by an exclusion list. Matching gene features and any other features
that overlap the deleted intervals are also removed. Coordinates of all
remaining features are adjusted accordingly.

Usage:
    python remove_transposases_and_phage.py input.gbff output.gbff
    python remove_transposases_and_phage.py input.gbff output.gbff --report removed.tsv
    python remove_transposases_and_phage.py input.gbff output.gbff --dry-run
    python remove_transposases_and_phage.py input.gbff output.gbff \
        --extra-pattern "relaxase" --extra-exclude "some product to keep"
"""
import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------------------------------
# Product matching
# --------------------------------------------------------------------------

# Case-insensitive regexes; a CDS is targeted if ANY of these match its product.
INCLUDE_PATTERNS = [
    # Transposases / mobile elements
    r"transposase",
    r"\bIS\d+\b",                      # IS1, IS3, IS200 ... family names
    r"insertion sequence",
    r"mobile element protein",
    r"transposon",
    r"\bTn\d+\b",                      # Tn3, Tn7, Tn916 ...
    r"\bTns[A-E]\b",                   # Tn7 machinery
    # Integrases / recombinases / resolvases (MGE-associated)
    r"integrase",
    r"recombinase",                    # tyrosine/serine recombinase, site-specific recombinase
    r"resolvase",                      # Tn3-family resolvases, serine resolvases
    r"invertase",                      # Hin/Gin/Min DNA invertases (serine recombinase family)
]

# Products that mention "phage" but are NOT phage-derived proteins
# (bacterial stress response, anti-phage defence systems, receptors, ...).
EXCLUDE_PATTERNS = [
    r"phage shock",                    # Psp system (PspA-G) - bacterial
    r"anti-?phage",
    r"phage defen[cs]e",
    r"phage exclusion",
    r"phage resistance",
    r"phage growth limitation",
    r"phage receptor",
    # Housekeeping recombinases that the broad 'recombinase' pattern would
    # otherwise catch. These are core chromosomal genes, not mobile elements.
    r"\brecA\b",                       # "recombinase RecA", "Protein RecA"
    r"recombinase A\b",
    r"\bXer[CD]\b",                    # chromosome dimer resolution (XerC/XerD)
    r"\bRecR\b|\bRecO\b|\bRecF\b",     # recombination mediator proteins
    r"integration host factor",        # IHF: not an integrase
]

_INCLUDE_RE = None
_EXCLUDE_RE = None


def compile_patterns(extra_include=None, extra_exclude=None):
    global _INCLUDE_RE, _EXCLUDE_RE
    inc = INCLUDE_PATTERNS + list(extra_include or [])
    exc = EXCLUDE_PATTERNS + list(extra_exclude or [])
    _INCLUDE_RE = re.compile("|".join(f"(?:{p})" for p in inc), re.IGNORECASE)
    _EXCLUDE_RE = re.compile("|".join(f"(?:{p})" for p in exc), re.IGNORECASE)


def is_target_product(product):
    """True if the product string should be excised."""
    if not product:
        return False
    if _EXCLUDE_RE.search(product):
        return False
    return bool(_INCLUDE_RE.search(product))


# --------------------------------------------------------------------------
# Interval / sequence helpers
# --------------------------------------------------------------------------

def get_feature_interval(feature):
    """Return (start, end) in 1-based inclusive coordinates for a feature."""
    loc = feature.location
    return int(loc.start) + 1, int(loc.end)


def merge_intervals(intervals):
    """Merge overlapping or adjacent intervals. Returns sorted list."""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(sorted_intervals[0])]
    for start, end in sorted_intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(m) for m in merged]


def sequence_without_intervals(seq, intervals):
    """Return new sequence string with given 1-based inclusive intervals removed."""
    if not intervals:
        return str(seq)
    kept = []
    prev_end = 0  # 0-based exclusive index
    for start, end in intervals:
        kept.append(str(seq[prev_end:start - 1]))
        prev_end = end
    kept.append(str(seq[prev_end:]))
    return "".join(kept)


def feature_overlaps_intervals(feature, intervals):
    """Check if a feature overlaps any of the given 1-based intervals."""
    loc = feature.location
    parts = loc.parts if isinstance(loc, CompoundLocation) else [loc]
    for part in parts:
        p_start = int(part.start) + 1
        p_end = int(part.end)
        for d_start, d_end in intervals:
            if not (p_end < d_start or p_start > d_end):
                return True
    return False


def adjust_location(loc, intervals):
    """Shift a location by subtracting lengths of intervals that end before it."""
    if isinstance(loc, CompoundLocation):
        new_parts = []
        for part in loc.parts:
            new_part = adjust_location(part, intervals)
            if new_part is None:
                return None
            new_parts.append(new_part)
        if len(new_parts) == 1:
            return new_parts[0]
        return CompoundLocation(new_parts, operator=loc.operator)

    f_start = int(loc.start) + 1
    shift = sum(d_end - d_start + 1 for d_start, d_end in intervals if d_end < f_start)
    return FeatureLocation(int(loc.start) - shift, int(loc.end) - shift, strand=loc.strand)


def process_record(record, deleted_intervals):
    """Create a new SeqRecord with specified intervals excised from the sequence."""
    deleted_intervals = merge_intervals(deleted_intervals)
    new_seq = Seq(sequence_without_intervals(record.seq, deleted_intervals))

    new_record = SeqRecord(
        new_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        annotations=record.annotations.copy(),
        dbxrefs=record.dbxrefs[:],
    )

    for feature in record.features:
        if feature.type == "source":
            new_loc = FeatureLocation(0, len(new_seq), strand=feature.location.strand)
            new_record.features.append(
                SeqFeature(location=new_loc, type="source",
                           qualifiers=feature.qualifiers.copy())
            )
            continue

        # Remove any feature that overlaps a deleted interval
        if feature_overlaps_intervals(feature, deleted_intervals):
            continue

        new_loc = adjust_location(feature.location, deleted_intervals)
        if new_loc is None:
            continue

        new_record.features.append(
            SeqFeature(location=new_loc, type=feature.type,
                       qualifiers=feature.qualifiers.copy())
        )

    return new_record


# --------------------------------------------------------------------------
# Main workflow
# --------------------------------------------------------------------------

def remove_target_sequence(input_gbff, output_gbff, report_path=None, dry_run=False):
    stats = {
        "records": 0,
        "records_with_deletions": 0,
        "target_cds": 0,
        "features_removed": 0,
        "features_kept": 0,
        "bp_removed": 0,
    }
    report_rows = []

    out_handle = None if dry_run else open(output_gbff, "w")
    try:
        for record in SeqIO.parse(input_gbff, "genbank"):
            target_intervals = []
            matching_genes = set()

            # Pass 1: identify target CDS intervals and their locus_tags
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                product = feature.qualifiers.get("product", [""])[0]
                if not is_target_product(product):
                    continue
                start, end = get_feature_interval(feature)
                target_intervals.append((start, end))
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                if locus_tag:
                    matching_genes.add(locus_tag)
                stats["target_cds"] += 1
                report_rows.append((record.id, locus_tag or "", start, end,
                                    feature.location.strand, product))

            # Pass 2: identify matching gene features by locus_tag or exact coordinates
            gene_intervals = []
            for feature in record.features:
                if feature.type != "gene":
                    continue
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                g_start, g_end = get_feature_interval(feature)
                if locus_tag in matching_genes:
                    gene_intervals.append((g_start, g_end))
                    continue
                for t_start, t_end in target_intervals:
                    if g_start == t_start and g_end == t_end:
                        gene_intervals.append((g_start, g_end))
                        break

            all_deleted = merge_intervals(target_intervals + gene_intervals)
            stats["records"] += 1

            if not all_deleted:
                if out_handle:
                    SeqIO.write(record, out_handle, "genbank")
                stats["features_kept"] += len(record.features)
                continue

            stats["records_with_deletions"] += 1
            stats["bp_removed"] += sum(end - start + 1 for start, end in all_deleted)

            new_record = process_record(record, all_deleted)
            if out_handle:
                SeqIO.write(new_record, out_handle, "genbank")

            stats["features_removed"] += len(record.features) - len(new_record.features)
            stats["features_kept"] += len(new_record.features)
    finally:
        if out_handle:
            out_handle.close()

    if report_path:
        with open(report_path, "w") as rep:
            rep.write("record\tlocus_tag\tstart\tend\tstrand\tproduct\n")
            for row in report_rows:
                rep.write("\t".join(str(x) for x in row) + "\n")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Excise transposase and phage CDS (and overlapping features) from a GenBank file."
    )
    parser.add_argument("input", help="input GenBank file (.gbff/.gbk)")
    parser.add_argument("output", help="output GenBank file")
    parser.add_argument("--report", help="write a TSV listing every excised CDS")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would be removed without writing the output GenBank")
    parser.add_argument("--extra-pattern", action="append", default=[],
                        help="additional case-insensitive regex to match in /product (repeatable)")
    parser.add_argument("--extra-exclude", action="append", default=[],
                        help="additional regex; products matching this are never removed (repeatable)")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    compile_patterns(args.extra_pattern, args.extra_exclude)

    print(f"Parsing {args.input} ...")
    stats = remove_target_sequence(args.input, args.output,
                                   report_path=args.report, dry_run=args.dry_run)

    if args.dry_run:
        print("Dry run: no output GenBank written.")
    else:
        print(f"Done. Wrote {stats['records']} records to {args.output}")
    print(f"  Records with deletions: {stats['records_with_deletions']}")
    print(f"  Transposase/phage CDS excised: {stats['target_cds']}")
    print(f"  Features removed (incl. genes/overlaps): {stats['features_removed']}")
    print(f"  Features kept: {stats['features_kept']}")
    print(f"  Total bp removed: {stats['bp_removed']}")
    if args.report:
        print(f"  Report written to {args.report}")


if __name__ == "__main__":
    main()