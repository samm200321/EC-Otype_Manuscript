#!/usr/bin/env python3
"""Filter multi-record GenBank files down to only the records (contigs) that
contain 3 or more ORFs each.

Unlike a whole-file filter, this looks at each record independently: a file
with 5 contigs where only 2 of them have >=3 CDS features will be rewritten
to contain just those 2 records. A record below the threshold is dropped
entirely (not truncated).

"ORF" here means a CDS feature; use --feature-type to count something else
(e.g. "gene").

Usage:
    # One output file per input file, records filtered, written to out_dir/
    python filter_genbank_records_by_orf_count.py input_dir/ out_dir/

    # Custom threshold
    python filter_genbank_records_by_orf_count.py input_dir/ out_dir/ --min-orfs 5

    # Individual files instead of a directory
    python filter_genbank_records_by_orf_count.py sample.gbff out_dir/

    # Drop whole files that end up with zero qualifying records
    python filter_genbank_records_by_orf_count.py input_dir/ out_dir/ --skip-empty-output

    # See per-record counts without writing anything
    python filter_genbank_records_by_orf_count.py input_dir/ out_dir/ --dry-run --report counts.tsv
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO

GBFF_EXTENSIONS = {".gbff", ".gbk", ".gb", ".genbank", ".gbf"}


def collect_input_files(inputs):
    """Expand a mix of directories and individual file paths into a file list."""
    files = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            for ext in GBFF_EXTENSIONS:
                files.extend(sorted(p.glob(f"*{ext}")))
        elif p.is_file():
            files.append(p)
        else:
            print(f"Warning: skipping missing path: {item}", file=sys.stderr)
    seen = set()
    unique_files = []
    for f in files:
        rp = f.resolve()
        if rp not in seen:
            seen.add(rp)
            unique_files.append(f)
    return unique_files


def filter_records(input_path, min_orfs, feature_type="CDS"):
    """Return (kept_records, per_record_counts) for one GenBank file.

    per_record_counts is a list of (record_id, orf_count, kept_bool).
    """
    kept = []
    counts = []
    for record in SeqIO.parse(input_path, "genbank"):
        n_orfs = sum(1 for f in record.features if f.type == feature_type)
        is_kept = n_orfs >= min_orfs
        counts.append((record.id, n_orfs, is_kept))
        if is_kept:
            kept.append(record)
    return kept, counts


def main():
    parser = argparse.ArgumentParser(
        description="Retain only GenBank records with >= N ORFs (CDS features) each, "
                    "from single- or multi-record GenBank files."
    )
    parser.add_argument(
        "inputs", nargs="+",
        help="one or more GenBank files and/or directories containing them "
             "(the last positional argument is treated as the output directory "
             "unless --output-dir is given explicitly)"
    )
    parser.add_argument("--output-dir",
                        help="directory to write filtered files into (overrides using "
                             "the last positional arg for this)")
    parser.add_argument("--min-orfs", type=int, default=3,
                        help="minimum number of ORFs a record must have to be kept (default: 3)")
    parser.add_argument("--feature-type", default="CDS",
                        help="feature type counted as an ORF (default: CDS)")
    parser.add_argument("--suffix", default="",
                        help="suffix to insert before the extension of output filenames, "
                             "e.g. --suffix .filtered -> sample.filtered.gbff (default: none, "
                             "same filename as input)")
    parser.add_argument("--skip-empty-output", action="store_true",
                        help="if a file ends up with zero qualifying records, don't write it at all")
    parser.add_argument("--report", help="write a TSV with one row per record: file, record id, "
                                          "ORF count, kept/dropped")
    parser.add_argument("--dry-run", action="store_true",
                        help="report per-record counts without writing any output files")
    args = parser.parse_args()

    if args.output_dir:
        output_dir = Path(args.output_dir)
        input_paths = args.inputs
    else:
        if len(args.inputs) < 2:
            parser.error("provide an output directory as the last argument, or use --output-dir")
        output_dir = Path(args.inputs[-1])
        input_paths = args.inputs[:-1]

    files = collect_input_files(input_paths)
    if not files:
        print("No GenBank files found in the given input(s).")
        sys.exit(1)

    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    report_rows = []
    total_records = 0
    total_kept_records = 0
    files_written = 0
    files_skipped_empty = 0
    files_errored = 0

    for f in files:
        try:
            kept_records, counts = filter_records(f, args.min_orfs, args.feature_type)
        except Exception as exc:
            print(f"Error: could not parse {f}: {exc}", file=sys.stderr)
            files_errored += 1
            report_rows.append((f.name, "ERROR", "", "", str(exc)))
            continue

        total_records += len(counts)
        total_kept_records += len(kept_records)
        for record_id, n_orfs, is_kept in counts:
            report_rows.append((f.name, record_id, n_orfs,
                                "kept" if is_kept else "dropped", ""))

        if not kept_records and args.skip_empty_output:
            files_skipped_empty += 1
            continue

        if not args.dry_run:
            if kept_records:
                stem = f.stem
                out_name = f"{stem}{args.suffix}{f.suffix}"
                out_path = output_dir / out_name
                SeqIO.write(kept_records, out_path, "genbank")
                files_written += 1
            elif not args.skip_empty_output:
                # write an empty-but-valid file only if user didn't ask to skip
                stem = f.stem
                out_name = f"{stem}{args.suffix}{f.suffix}"
                out_path = output_dir / out_name
                open(out_path, "w").close()
                files_written += 1
        else:
            if kept_records or not args.skip_empty_output:
                files_written += 1

    if args.report:
        with open(args.report, "w") as rep:
            rep.write(f"file\trecord_id\t{args.feature_type}_count\tstatus\terror\n")
            for row in report_rows:
                rep.write("\t".join(str(x) for x in row) + "\n")

    action = "Would write" if args.dry_run else "Wrote"
    print(f"Scanned {len(files)} file(s), {total_records} record(s) total, "
          f"threshold: >= {args.min_orfs} {args.feature_type} feature(s) per record")
    print(f"  Kept {total_kept_records}/{total_records} record(s)")
    print(f"  {action} {files_written} output file(s) to {output_dir}")
    if args.skip_empty_output:
        print(f"  Skipped {files_skipped_empty} file(s) with zero qualifying records")
    if files_errored:
        print(f"  Failed to parse {files_errored} file(s) (see stderr)")
    if args.report:
        print(f"  Report written to {args.report}")


if __name__ == "__main__":
    main()