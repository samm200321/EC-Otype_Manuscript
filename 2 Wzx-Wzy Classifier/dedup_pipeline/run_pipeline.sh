#!/usr/bin/env bash
# De-duplication pipeline:
#   1) exclude matched-labels queries from data/all_prot_2.faa
#   2) exclude proteins in data/short_proteins_lt100aa.csv
#   3) mmseqs2 all-vs-all search + identity histogram -> pick threshold
#   4) mmseqs2 cluster at MIN_SEQ_ID -> final de-duplicated fasta
#
# Usage: ./run_pipeline.sh [MIN_SEQ_ID]   (default 0.9)
set -euo pipefail
cd "$(dirname "$0")"

MIN_SEQ_ID="${1:-0.9}"

python3 01_filter_exclusions.py
python3 02_pairwise_identities.py
python3 03_plot_identities.py
python3 04_cluster_dedup.py --min-seq-id "$MIN_SEQ_ID"

echo "Pipeline complete. Results in dedup_pipeline/results/"
