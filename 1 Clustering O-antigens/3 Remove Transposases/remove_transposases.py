#!/usr/bin/env python3
"""Parse a GenBank file and physically delete the nucleotide sequence of CDS
features where the /product contains 'transposase'. Matching gene features
and any other features that overlap the deleted intervals are also removed.
Coordinates of all remaining features are adjusted accordingly.

Usage:
    python remove_transposases.py input.gbff output.gbff
"""
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
        # Keep seq[prev_end : start-1]
        kept.append(str(seq[prev_end:start - 1]))
        prev_end = end  # end is 1-based inclusive, next kept starts here (0-based)
    kept.append(str(seq[prev_end:]))
    return "".join(kept)


def feature_overlaps_intervals(feature, intervals):
    """Check if a feature overlaps any of the given 1-based intervals."""
    loc = feature.location
    if isinstance(loc, CompoundLocation):
        for part in loc.parts:
            p_start = int(part.start) + 1
            p_end = int(part.end)
            for d_start, d_end in intervals:
                if not (p_end < d_start or p_start > d_end):
                    return True
        return False
    else:
        f_start = int(loc.start) + 1
        f_end = int(loc.end)
        for d_start, d_end in intervals:
            if not (f_end < d_start or f_start > d_end):
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
    f_end = int(loc.end)

    shift = sum(d_end - d_start + 1 for d_start, d_end in intervals if d_end < f_start)
    new_start = int(loc.start) - shift
    new_end = int(loc.end) - shift
    return FeatureLocation(new_start, new_end, strand=loc.strand)


def process_record(record, deleted_intervals):
    """Create a new SeqRecord with specified intervals excised from the sequence."""
    deleted_intervals = merge_intervals(deleted_intervals)
    new_seq_str = sequence_without_intervals(record.seq, deleted_intervals)
    new_seq = Seq(new_seq_str)

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
            new_feature = SeqFeature(
                location=new_loc,
                type="source",
                qualifiers=feature.qualifiers.copy(),
            )
            new_record.features.append(new_feature)
            continue

        # Remove any feature that overlaps a deleted interval
        if feature_overlaps_intervals(feature, deleted_intervals):
            continue

        new_loc = adjust_location(feature.location, deleted_intervals)
        if new_loc is None:
            continue

        new_feature = SeqFeature(
            location=new_loc,
            type=feature.type,
            qualifiers=feature.qualifiers.copy(),
        )
        new_record.features.append(new_feature)

    return new_record


def remove_transposase_sequence(input_gbff, output_gbff):
    records_processed = 0
    total_features_removed = 0
    total_features_kept = 0
    total_bp_removed = 0
    records_with_deletions = 0

    with open(output_gbff, "w") as out_handle:
        for record in SeqIO.parse(input_gbff, "genbank"):
            transposase_intervals = []
            matching_genes = set()

            # Pass 1: identify transposase CDS intervals and their locus_tags
            for feature in record.features:
                if feature.type == "CDS":
                    product = feature.qualifiers.get("product", [""])[0]
                    if "transposase" in product.lower():
                        start, end = get_feature_interval(feature)
                        transposase_intervals.append((start, end))
                        locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                        if locus_tag:
                            matching_genes.add(locus_tag)

            # Pass 2: identify matching gene features by locus_tag or exact coordinates
            gene_intervals = []
            for feature in record.features:
                if feature.type == "gene":
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    g_start, g_end = get_feature_interval(feature)
                    if locus_tag in matching_genes:
                        gene_intervals.append((g_start, g_end))
                        continue
                    for t_start, t_end in transposase_intervals:
                        if g_start == t_start and g_end == t_end:
                            gene_intervals.append((g_start, g_end))
                            if locus_tag:
                                matching_genes.add(locus_tag)
                            break

            all_deleted = merge_intervals(transposase_intervals + gene_intervals)

            if not all_deleted:
                SeqIO.write(record, out_handle, "genbank")
                records_processed += 1
                total_features_kept += len(record.features)
                continue

            records_with_deletions += 1
            bp_removed = sum(end - start + 1 for start, end in all_deleted)
            total_bp_removed += bp_removed

            new_record = process_record(record, all_deleted)
            SeqIO.write(new_record, out_handle, "genbank")
            records_processed += 1

            removed = len(record.features) - len(new_record.features)
            total_features_removed += removed
            total_features_kept += len(new_record.features)

    return {
        "records": records_processed,
        "records_with_deletions": records_with_deletions,
        "features_removed": total_features_removed,
        "features_kept": total_features_kept,
        "bp_removed": total_bp_removed,
    }


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.gbff output.gbff")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if not Path(input_path).exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    print(f"Parsing {input_path} ...")
    stats = remove_transposase_sequence(input_path, output_path)

    print(f"Done. Wrote {stats['records']} records to {output_path}")
    print(f"  Records with deletions: {stats['records_with_deletions']}")
    print(f"  Features removed: {stats['features_removed']}")
    print(f"  Features kept: {stats['features_kept']}")
    print(f"  Total bp removed: {stats['bp_removed']}")


if __name__ == "__main__":
    main()
