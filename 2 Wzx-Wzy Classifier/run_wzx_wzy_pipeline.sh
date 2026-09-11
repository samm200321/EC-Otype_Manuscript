#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Assign O-antigen locus genes to wzx/wzy by homology to a labeled reference
# set (>=95% identity), then de-duplicate the unmatched remainder by
# clustering it on its own.
#
# Usage:
#   ./run_wzx_wzy_pipeline.sh <labeled.fasta> <unlabeled.fasta> <outdir>
#
# Requirements: mmseqs2 (https://github.com/soedinglab/MMseqs2), python3+pandas
# ---------------------------------------------------------------------------
set -euo pipefail

LABELED_FASTA=${1:?Usage: $0 <labeled.fasta> <unlabeled.fasta> <outdir>}
UNLABELED_FASTA=${2:?Usage: $0 <labeled.fasta> <unlabeled.fasta> <outdir>}
OUTDIR=${3:?Usage: $0 <labeled.fasta> <unlabeled.fasta> <outdir>}

THREADS=${THREADS:-8}
ID_THRESHOLD=${ID_THRESHOLD:-0.85}      # fraction, i.e. 0.95 = 95% identity
COV_THRESHOLD=${COV_THRESHOLD:-0.99}     # min alignment coverage (of the shorter/target seq)
CLUSTER_ID=${CLUSTER_ID:-0.90}          # identity threshold for de-duplicating the "other" set

mkdir -p "$OUTDIR"
TMP="$OUTDIR/tmp"
mkdir -p "$TMP"

echo "== 1. Building mmseqs databases =="
mmseqs createdb "$LABELED_FASTA"   "$OUTDIR/labeledDB"
mmseqs createdb "$UNLABELED_FASTA" "$OUTDIR/unlabeledDB"

echo "== 2. Searching unlabeled genes against labeled wzx/wzy set =="
# --max-seqs kept modest (10) so we can pick the true best-identity hit per
# query in postprocessing, rather than trusting search's default e-value sort.
mmseqs search "$OUTDIR/unlabeledDB" "$OUTDIR/labeledDB" "$OUTDIR/resultDB" "$TMP" \
    --min-seq-id "$ID_THRESHOLD" \
    -c "$COV_THRESHOLD" --cov-mode 0 \
    --max-seqs 10 \
    -s 7.5 \
    --threads "$THREADS"

echo "== 3. Converting search results to TSV =="
mmseqs convertalis "$OUTDIR/unlabeledDB" "$OUTDIR/labeledDB" "$OUTDIR/resultDB" "$OUTDIR/hits.tsv" \
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits" \
    --threads "$THREADS"

echo "== 4. Extracting all labeled and unlabeled sequence IDs =="
mmseqs createtsv "$OUTDIR/labeledDB"   "$OUTDIR/labeledDB"   "$OUTDIR/labeled_ids.tsv"
mmseqs createtsv "$OUTDIR/unlabeledDB" "$OUTDIR/unlabeledDB" "$OUTDIR/unlabeled_ids.tsv"

echo "== 5. Postprocessing: assign labels from best hits (python) =="
python3 "$(dirname "$0")/assign_labels.py" \
    --hits "$OUTDIR/hits.tsv" \
    --labeled-fasta "$LABELED_FASTA" \
    --unlabeled-fasta "$UNLABELED_FASTA" \
    --id-threshold "$ID_THRESHOLD" \
    --cov-threshold "$COV_THRESHOLD" \
    --outdir "$OUTDIR"

echo "== 6. Clustering the unmatched ('other') sequences to de-duplicate =="
if [ -s "$OUTDIR/unmatched.fasta" ]; then
    mmseqs easy-cluster "$OUTDIR/unmatched.fasta" "$OUTDIR/other_clusters" "$TMP" \
        --min-seq-id "$CLUSTER_ID" \
        -c 0.8 --cov-mode 0 \
        --threads "$THREADS"
    # easy-cluster produces: other_clusters_cluster.tsv (rep<TAB>member), other_clusters_rep_seq.fasta
else
    echo "No unmatched sequences to cluster - skipping."
fi

echo "== 7. Building final assignment table =="
python3 "$(dirname "$0")/build_final_table.py" \
    --unlabeled-fasta "$UNLABELED_FASTA" \
    --matched "$OUTDIR/matched_labels.tsv" \
    --cluster-tsv "$OUTDIR/other_clusters_cluster.tsv" \
    --out "$OUTDIR/final_assignments.tsv"

echo "Done. See $OUTDIR/final_assignments.tsv"
