# Create MMseqs2 database
mmseqs createdb oag_prot.faa proteinDB

# Cluster (linclust is fast; cluster is more sensitive)
mmseqs cluster proteinDB clusterDB tmp \
    --min-seq-id 0.8 \
    --cov-mode 0 \
    -c 0.8 \
    --cluster-mode 2 \
    -e 0.001 \
    --threads 8

# Extract cluster TSV (representative → member)
mmseqs createtsv proteinDB proteinDB clusterDB clusters.tsv

# Get representative sequences (the "pan-proteome")
mmseqs result2repseq proteinDB clusterDB repDB
mmseqs convert2fasta repDB representative_seqs.faa