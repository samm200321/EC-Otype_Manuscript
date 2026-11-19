#!/usr/bin/env python3
"""Embed every protein in a multi-FASTA file with ESM-C 600M (ESM++), every layer.

For each protein, one forward pass is run with output_hidden_states=True and the
hidden state of every layer is mean-pooled over real residues (excluding
BOS/EOS/padding): state 0 is the token-embedding output, states 1..n-1 are the
outputs of each transformer block, and the final state is the post-LayerNorm
last_hidden_state. For Synthyra/ESMplusplus_large (ESM-C 600M-equivalent) that
is 36 + 1 = 37 states of dimension 1152 per protein.

Output is a wide CSV: one row per sequence, columns `layer{L}_dim_{D}`.

Usage:
    python esmc600m_embeddings.py proteins.faa -o embeddings.csv
"""

import argparse
import os
import sys

# Must be set before torch import (MPS op coverage is not 100% complete on macOS;
# unimplemented ops silently run on CPU instead of raising).
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from tqdm.auto import tqdm
from transformers import AutoConfig, AutoModelForMaskedLM

DEFAULT_MODEL_NAME = "Synthyra/ESMplusplus_large"  # ESM-C 600M-equivalent, transformers-compatible
DEFAULT_MAX_LENGTH = 2048  # ESM-C context window limit (tokens); longer sequences are truncated


def parse_args():
    p = argparse.ArgumentParser(
        description="Mean-pool ESM-C 600M (ESM++) embeddings for every layer of each protein "
                    "in a multi-FASTA file and save them to a CSV."
    )
    p.add_argument("fasta", help="Input multi-FASTA protein file")
    p.add_argument(
        "-o", "--out",
        default=None,
        help="Output CSV path (default: <fasta_basename>_all_layers_embeddings.csv)",
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
        help="Sequences per embedding batch; start small on MPS, increase if memory allows (default: 2)",
    )
    p.add_argument(
        "--device",
        choices=["auto", "cpu", "mps", "cuda"],
        default="auto",
        help="Compute device (default: auto -> cuda > mps > cpu)",
    )
    return p.parse_args()


def pick_device(choice):
    if choice != "auto":
        return torch.device(choice)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_fasta(path, max_length):
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).upper().replace("*", "")
        if not seq:
            print(f"Warning: skipping {rec.id} (empty sequence)", file=sys.stderr)
            continue
        if len(seq) > max_length:
            print(f"Warning: truncating {rec.id} from {len(seq)} to {max_length} aa", file=sys.stderr)
            seq = seq[:max_length]
        records.append((rec.id, seq))
    return records


def load_embedder(model_name, device):
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.attn_backend = "sdpa"

    model = AutoModelForMaskedLM.from_pretrained(
        model_name,
        config=config,
        trust_remote_code=True,
    ).to(device).eval()

    return model, model.tokenizer


def mean_pool_hidden_states(hidden_states, attention_mask):
    """Mean-pool per-residue embeddings, excluding padding and BOS/EOS tokens.

    hidden_states: (batch, seq_len, hidden_dim)
    attention_mask: (batch, seq_len) - 1 for real tokens, 0 for padding
    """
    mask = attention_mask.clone()
    mask[:, 0] = 0  # BOS
    seq_lengths = attention_mask.sum(dim=1)  # includes BOS+EOS+residues
    for i, length in enumerate(seq_lengths):
        eos_pos = int(length.item()) - 1
        if eos_pos >= 0:
            mask[i, eos_pos] = 0  # EOS

    mask = mask.unsqueeze(-1).to(hidden_states.dtype)  # (batch, seq_len, 1)
    summed = (hidden_states * mask).sum(dim=1)
    counts = mask.sum(dim=1).clamp(min=1)
    return (summed / counts).cpu().float().numpy()


def embed_all_layers(records, model, tokenizer, device, max_length, batch_size):
    """Mean-pool every layer's hidden states for all records.

    Returns (ids, per_layer_embeddings, failed_ids) where per_layer_embeddings
    is a list (one entry per layer state) of 2-D arrays (n_seqs, hidden_dim).
    """
    all_ids, failed = [], []
    per_layer_arrays = None  # sized on the first successful batch
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

            hidden_states = output.hidden_states  # tuple of (batch, seq_len, hidden_dim)
            if per_layer_arrays is None:
                per_layer_arrays = [[] for _ in range(len(hidden_states))]

            for layer_idx, layer_hidden in enumerate(hidden_states):
                pooled = mean_pool_hidden_states(layer_hidden, inputs["attention_mask"])
                per_layer_arrays[layer_idx].append(pooled)

            all_ids.extend(batch_ids)
        except Exception as e:
            print(f"Warning: failed to embed batch {batch_ids}: {e}", file=sys.stderr)
            failed.extend(batch_ids)

    per_layer_embeddings = [
        np.vstack(layer) if layer else np.empty((0, 0), dtype=np.float32)
        for layer in (per_layer_arrays or [])
    ]
    return all_ids, per_layer_embeddings, failed


def save_embeddings_csv(ids, per_layer_embeddings, out_path):
    num_states = len(per_layer_embeddings)
    embedding_dim = per_layer_embeddings[0].shape[1]

    frames = [pd.DataFrame({"id": ids})]
    for layer_idx in range(num_states):
        col_names = [f"layer{layer_idx}_dim_{i}" for i in range(embedding_dim)]
        frames.append(pd.DataFrame(per_layer_embeddings[layer_idx], columns=col_names))

    df = pd.concat(frames, axis=1)
    df.to_csv(out_path, index=False)
    return df


def main():
    args = parse_args()

    device = pick_device(args.device)
    print(f"Using device: {device}")

    records = load_fasta(args.fasta, args.max_length)
    if not records:
        sys.exit(f"No usable protein sequences found in {args.fasta}")
    print(f"Loaded {len(records)} sequences from {args.fasta}")

    print(f"Loading embedding model {args.model_name} ...")
    model, tokenizer = load_embedder(args.model_name, device)

    ids, per_layer_embeddings, failed = embed_all_layers(
        records, model, tokenizer, device, args.max_length, args.batch_size
    )
    if failed:
        print(f"Warning: {len(failed)} sequences failed to embed: {failed}", file=sys.stderr)
    if not ids:
        sys.exit("No sequences could be embedded; nothing to save.")

    num_states = len(per_layer_embeddings)
    embedding_dim = per_layer_embeddings[0].shape[1]
    print(f"Embedded {len(ids)} / {len(records)} sequences across {num_states} layer states "
          f"({embedding_dim} dims each)")

    out = args.out or os.path.splitext(args.fasta)[0] + "_all_layers_embeddings.csv"
    df = save_embeddings_csv(ids, per_layer_embeddings, out)
    print(f"Saved {df.shape[0]} sequences x {num_states} layers "
          f"({df.shape[1] - 1} feature columns total) to {out}")


if __name__ == "__main__":
    main()
