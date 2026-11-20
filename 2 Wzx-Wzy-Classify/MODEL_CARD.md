# Model card — `wzx-wzy-rf.joblib`

Random forest classifier that labels each protein as **Wzx** (O-antigen flippase), **Wzy**
(O-antigen polymerase) or **Other**, from ESM-C protein embeddings.

## Model details

- **Type:** `StandardScaler` + `RandomForestClassifier` (200 trees, `random_state=42`), bundled with its fitted `LabelEncoder` as a joblib dict: `{"pipeline": ..., "label_encoder": ...}`.
- **File:** `2 Wzx-Wzy-Classify/wzx-wzy-rf.joblib` (~360 KB), produced by `ESMC_training.ipynb` in October 2026.
- **Features:** the 1152-dim mean-pooled hidden state of layer 11 (transformer block 11) of ESM-C 600M (`Synthyra/ESMplusplus_large`), pooled over real residues only (BOS/EOS and padding excluded). Sequences longer than 2,048 tokens are truncated.

## Training data

- Source: `oag_wzx_wzy_dedup.faa` — 1,722 deduplicated proteins from curated *E. coli* O-antigen loci.
- Labels are derived from the protein IDs: IDs containing `wzx` → Wzx, `wzy` → Wzy, everything else → Other.
- Removed before training: proteins with `novel` IDs, the 25 short proteins in `data/short_proteins.csv`, and the 369 manually excluded proteins in `data/exclusion list.csv`.
- Final set: **1,634 proteins** — 166 Wzx, 165 Wzy, 1,303 Other.
- Split: 80/20 stratified (`random_state=40`) → 1,307 training / 327 held-out test.

## Training procedure

- All 37 ESM-C layer states (0–36) were evaluated per model with 5-fold stratified cross-validation (macro-F1) on the training split; **layer 11** gave the best random forest result (CV macro-F1 0.997).
- KNN and SVM performed comparably (macro-F1 ≈ 0.99); the random forest was selected and exported.
- Embeddings were generated with `esmc600m_embeddings.py`; the exported pipeline includes the fitted scaler so inference preprocessing matches training exactly.

## Evaluation (held-out test split, n = 327)

| Class  | Support | Precision | Recall | F1   |
|--------|---------|-----------|--------|------|
| Other  | 261     | 1.00      | 1.00   | 1.00 |
| Wzx    | 33      | 1.00      | 0.97   | 0.98 |
| Wzy    | 33      | 0.97      | 1.00   | 0.99 |

- Accuracy: 326/327 (99.7%); macro-F1: 0.99. The single error is one Wzx protein predicted as Wzy.
- The confusion matrix is rendered in `ESMC_training.ipynb`.

## Intended use

Classifying wzx/wzy genes within *E. coli* O-antigen (rfb/wb*) loci, e.g. for database curation
or genotyping workflows.

## How to use

```bash
python classify_wzx_wzy.py proteins.faa -o predictions.csv
```

The script embeds every input protein with ESM-C 600M layer 11 and writes one row per protein
with the predicted class, confidence and per-class probabilities. The ESM-C model is
downloaded from HuggingFace on first use. A different embedding model or layer is rejected by
the script's feature-dimension check.

## Limitations

- **Other is a catch-all class**, dominated by O-antigen biosynthesis proteins; predicted
  probabilities are meaningful only for proteins from this context.
- **Class imbalance:** Wzx and Wzy each make up ~10% of the training set; metrics are
  macro-averaged for this reason.
- **Labels are inherited from database annotations**, not per-protein experimental validation.
- **Embedding-model dependency:** predictions require the same `Synthyra/ESMplusplus_large`
  weights; changes to the upstream model may shift predictions.
