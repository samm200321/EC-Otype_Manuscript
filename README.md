# GEOD Manuscript

This repo contains the python notebooks and scripts required to reproduce the key results and figures in *A gene cluster-based O-antigen typing database for Escherichia coli*

**Folder 1 - Clustering of O-antigens**
Each sub-directory should be run sequentially. 
This workflow implements the screening, pre-processing and clustering. 
This produced 2 GenBank files, containing established O-antigen loci and putative novel loci. Both are manually curated.


**Folder 2 - Wzx-Wzy Classifier**
This folder implements the Wzx-Wzy classifier. The fasta file of all proteins used for training is provided. To reproduce, embeddings would need to be regenerated using the esmc600m_embeddings.py. See MODEL_CARD.md for details of the trained model (`wzx-wzy-rf.joblib`).


**Folder 3 - O-antigen Pan-Genome**
This folder reproduces the statistics and figures about the O-antigen pangenome.

**Folder 4 - Benchmarking**
This contains the raw output of Kaptive, EcOH and ECTyper. The notebook performs the statistical comparisons.

## How to run

Notebooks use relative paths, so run each one from inside its own folder. Large GenBank files (`*.gbff`) are stored with Git LFS — run `git lfs install` before cloning.

**Prerequisites:** Python 3 with pandas, Biopython, scikit-learn, seaborn and plotly, plus BLAST+ (folder 1), PPanGGOLiN (sub-folder 1.4), MMseqs2 (folder 3) and PyTorch/transformers (folder 2).

**Folder 1** — run the sub-folders in order:
1. *1 Extraction of galf-gnd loci*: `extract_oag_loci.py` locates each O-antigen locus between galF and gnd using BLAST+. Requires E. coli assemblies in an `ec_refseq/` folder (not included in the repo).
2. *2 Remove IS Elements*: `process_isfinder.ipynb` excises the ISfinder BLAST hits from the loci.
3. *3 Remove Transposases*: `remove_transposases_and_phage.py` (or `remove_transposases.py`) excises transposase/phage CDS; `filter_min_orfs.py` keeps loci with >= 3 ORFs.
4. *4 Build Pangenome Matrix*: `filter_genbank_by_orf_count.py`, `reverse_complement_genbank.py` and `split_gbff_for_ppanggolin.py` prepare the per-locus GenBank files and the manifest `ppanggolin_input.tsv`; run PPanGGOLiN on the manifest to build the presence/absence matrix (`matrix.csv`).
5. *5 O-antigen Loci Clustering*: `process_clusters.ipynb`, then `extract_established.ipynb` and `novel_oag_validation.ipynb`. These read the cleaned loci file `O-antigen_db.transposase_clean.min3orfs_rc.gbff` produced by steps 1-4 (not included in the repo).
6. *6 Process Oagdb*: `process_oag_db.ipynb` combines the loci into the final database.

**Folder 2** — regenerate the embeddings, then run `ESMC_training.ipynb` (it loads the embeddings from `data/`; the two curation lists it also reads are not included):

```bash
python esmc600m_embeddings.py oag_wzx_wzy_dedup.faa
```

The trained model is included, so new proteins can be classified directly:

```bash
python classify_wzx_wzy.py proteins.faa -o predictions.csv
```

The ESM-C model (Synthyra/ESMplusplus_large) is downloaded from HuggingFace on first use.

**Folder 3** — run `process_pan.ipynb` (uses the MMseqs2 clusters in `clustering/clusters.tsv`; the clustering command is in `clustering/cluster.sh`) and `Gene_stats.ipynb` (uses `EC_Oag.gbk`). `pairwise_id.py` computes all-vs-all protein identities with MMseqs2.

**Folder 4** — run `mcnemar_geod_benchmark.ipynb` on `benchmark_summary.csv`. The raw Kaptive, EcOH and ECTyper outputs are included for reference.
