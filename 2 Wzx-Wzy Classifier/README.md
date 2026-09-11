## To run this workflow follow the following steps

1) Use the esmc600m_mac_embeddings.ipynb to generate embeddings for the proteins in data/all_prot_2.faa, wzx_proteins.faa, wzx_proteins.faa
2) Run the following command to identify Wzx/Wzy in all_prot_2.faa:
```bash 
run_wzx_wzy_pipeline.sh data/wzx_wzy.faa all_prot_2.faa
```
3) Run the build_exclusion_list.py
4) Run the dedup_pipeline in dedup_pipeline/run_pipeline.sh
5) Build the classifier in the Wzx-Wzy-MLP.ipynb notebook