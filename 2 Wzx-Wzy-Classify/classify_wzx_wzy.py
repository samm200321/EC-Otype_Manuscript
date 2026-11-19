#!/usr/bin/env python3
"""Classify every protein in a multi-FASTA file as "Wzx", "Wzy", or "Other".

Each protein is embedded with ESM-C 600M (Synthyra/ESMplusplus_large) and the
layer-11 hidden state is mean-pooled over real residues (excluding BOS/EOS and
padding), exactly as in esmc600m_embeddings.py. The 1152-dim embedding is then
classified with the pre-trained pipeline in wzx-wzy-rf.joblib, a dict bundling
{'pipeline': StandardScaler + RandomForestClassifier, 'label_encoder': LabelEncoder}.

Output is a CSV table with one row per protein: id, description, length,
prediction (Wzx/Wzy/Other), the prediction confidence and the per-class
probabilities from the forest.

Usage:
    python classify_wzx_wzy.py proteins.faa -o predictions.csv
"""

import argparse
import os
import sys

import joblib
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from tqdm.auto import tqdm

# Reuse the exact embedding pipeline that generated the training embeddings so
# preprocessing is identical (same cleaning, truncation, pooling, layer indexing).
from esmc600m_embeddings import (
    DEFAULT_MODEL_NAME,
    DEFAULT_MAX_LENGTH,
    load_embedder,
    load_fasta,
    mean_pool_hidden_states,
    pick_device,
)

# Layer-11 hidden state the classifier was trained on (state 0 = token-embedding
# output, states 1..n = transformer blocks).
LAYER = 11


def parse_args():
    here = os.path.dirname(os.path.abspath(__file__))
    default_classifier = os.path.join(here, "wzx-wzy-rf.joblib")

    p = argparse.ArgumentParser(
        description="Predict Wzx/Wzy/Other for every protein in a multi-FASTA file "
                    "using layer-11 ESM-C 600M embeddings and a pre-trained random forest."
    )
    p.add_argument("fasta", help="Input multi-FASTA protein file")
    p.add_argument(
        "-o", "--out",
        default=None,
        help="Output table path (default: <fasta_basename>_wzx_wzy_predictions.csv)",
    )
    p.add_argument(
        "-c", "--classifier",
        default=default_classifier,
        help=f"Trained model joblib file: dict with 'pipeline' (scaler + RF) and "
             f"'label_encoder' (default: {default_classifier})",
    )
    p.add_argument(
        "--model-name",
        default=DEFAULT_MODEL_NAME,
        help=f"HuggingFace model to embed with (default: {DEFAULT_MODEL_NAME})",
    )
    p.add_argument(
        "--max-length",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help=f"Max tokenized length; longer sequences are truncated (default: {DEFAULT_MAX_LENGTH})",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=2,
        help="Sequences per embedding batch; increase if memory allows (default: 2)",
    )
    p.add_argument(
        "--device",
        choices=["auto", "cpu", "mps", "cuda"],
        default="auto",
        help="Compute device (default: auto -> cuda > mps > cpu)",
    )
    p.add_argument(
        "--no-probs",
        action="store_true",
        help="Omit per-class probability columns from the output table",
    )
    return p.parse_args()


def load_descriptions(path):
    """Map sequence id -> full FASTA header line, for the output table."""
    return {rec.id: rec.description for rec in SeqIO.parse(path, "fasta")}


def embed_layer(records, model, tokenizer, device, max_length, batch_size, layer):
    """Mean-pool one layer's hidden states for all records.

    Returns (ids, X, failed_ids) where X is (n_embedded, hidden_dim).
    """
    all_ids, failed = [], []
    pooled_batches = []
    batches = [records[i:i + batch_size] for i in range(0, len(records), batch_size)]

    for batch in tqdm(batches, desc=f"Embedding sequences ({device})"):
        batch_ids = [r[0] for r in batch]
        batch_seqs = [r[1] for r in batch]
        try:
            inputs = tokenizer(
                batch_seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_length,
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}

            with torch.inference_mode():
                output = model(**inputs, output_hidden_states=True)

            pooled = mean_pool_hidden_states(output.hidden_states[layer], inputs["attention_mask"])
            pooled_batches.append(pooled)
            all_ids.extend(batch_ids)
        except Exception as e:
            print(f"Warning: failed to embed batch {batch_ids}: {e}", file=sys.stderr)
            failed.extend(batch_ids)

    X = np.vstack(pooled_batches) if pooled_batches else np.empty((0, 0), dtype=np.float32)
    return all_ids, X, failed


def load_model(path):
    """Load the joblib bundle saved as {'pipeline': ..., 'label_encoder': ...}."""
    bundle = joblib.load(path)
    if not (isinstance(bundle, dict) and "pipeline" in bundle and "label_encoder" in bundle):
        sys.exit(f"Error: {path} must be a dict with 'pipeline' (scaler + classifier) "
                 "and 'label_encoder' keys.")
    pipeline = bundle["pipeline"]
    label_encoder = bundle["label_encoder"]
    if not hasattr(label_encoder, "classes_"):
        sys.exit(f"Error: 'label_encoder' in {path} is not a fitted sklearn LabelEncoder.")
    return pipeline, label_encoder


def feature_frame(X, ids, pipeline):
    """Wrap the pooled layer-11 embeddings in a DataFrame whose columns (names and
    order) match the feature names the scaler was fitted with."""
    feature_names = None
    scaler = pipeline.steps[0][1]
    if hasattr(scaler, "feature_names_in_"):
        feature_names = list(scaler.feature_names_in_)

    expected = getattr(pipeline, "n_features_in_", None)
    if feature_names is not None and len(feature_names) != X.shape[1]:
        sys.exit(f"Error: model was fitted with {len(feature_names)} features "
                 f"({feature_names[0]}...) but embeddings have {X.shape[1]}; "
                 f"check the embedding model/layer.")
    if feature_names is None:
        feature_names = [f"layer{LAYER}_dim_{i}" for i in range(X.shape[1])]

    if expected is not None and expected != X.shape[1]:
        sys.exit(f"Error: model expects {expected} features but got {X.shape[1]}; "
                 f"check the embedding model/layer.")
    return pd.DataFrame(X, columns=feature_names, index=ids)


def classify(X, pipeline, label_encoder):
    """Return (predicted class names, per-class probabilities, probability column names)."""
    proba = pipeline.predict_proba(X)
    predicted = pipeline.predict(X)
    pred_names = label_encoder.inverse_transform(predicted)
    proba_names = label_encoder.inverse_transform(pipeline.classes_)
    return list(pred_names), proba, list(proba_names)


def build_table(ids, records_meta, descriptions, pred_names, proba, class_names, with_probs):
    lengths = dict(records_meta)
    df = pd.DataFrame({
        "id": ids,
        "description": [descriptions.get(i, "") for i in ids],
        "length": [lengths.get(i, np.nan) for i in ids],
        "prediction": pred_names,
    })

    if with_probs and proba is not None:
        prob_cols = [f"prob_{name}" for name in class_names]
        df[prob_cols] = proba
        df["confidence"] = df[prob_cols].max(axis=1)

    return df


def main():
    args = parse_args()

    device = pick_device(args.device)
    print(f"Using device: {device}")

    records = load_fasta(args.fasta, args.max_length)
    if not records:
        sys.exit(f"No usable protein sequences found in {args.fasta}")
    descriptions = load_descriptions(args.fasta)
    print(f"Loaded {len(records)} sequences from {args.fasta}")

    print(f"Loading embedding model {args.model_name} ...")
    model, tokenizer = load_embedder(args.model_name, device)

    ids, X, failed = embed_layer(
        records, model, tokenizer, device, args.max_length, args.batch_size, layer=LAYER
    )
    if failed:
        print(f"Warning: {len(failed)} sequences failed to embed and get no prediction: {failed}",
              file=sys.stderr)
    if not ids:
        sys.exit("No sequences could be embedded; nothing to classify.")

    print(f"Embedded {len(ids)} / {len(records)} sequences ({X.shape[1]} features each)")

    print(f"Loading model {args.classifier} ...")
    pipeline, label_encoder = load_model(args.classifier)
    X = feature_frame(X, ids, pipeline)

    pred_names, proba, proba_names = classify(X, pipeline, label_encoder)
    df = build_table(ids, [(r[0], len(r[1])) for r in records], descriptions,
                     pred_names, proba, proba_names, with_probs=not args.no_probs)

    out = args.out or os.path.splitext(args.fasta)[0] + "_wzx_wzy_predictions.csv"
    df.to_csv(out, index=False)

    counts = pd.Series(pred_names).value_counts()
    summary = ", ".join(f"{name}: {count}" for name, count in counts.items())
    print(f"Wrote {len(df)} predictions to {out} ({summary})")


if __name__ == "__main__":
    main()
