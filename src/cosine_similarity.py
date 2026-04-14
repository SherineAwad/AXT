#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse

warnings.filterwarnings('ignore')

def calculate_cosine_similarity_sparse(X_ref, X_query, ref_labels, query_labels):
    ref_celltypes = np.unique(ref_labels)
    query_celltypes = np.unique(query_labels)
    similarity_matrix = np.zeros((len(ref_celltypes), len(query_celltypes)))

    ref_profiles = {}
    for ct in ref_celltypes:
        ref_mask = (ref_labels == ct)
        if np.sum(ref_mask) > 0:
            if sparse.issparse(X_ref):
                ref_profiles[ct] = np.array(X_ref[ref_mask].mean(axis=0)).flatten()
            else:
                ref_profiles[ct] = X_ref[ref_mask].mean(axis=0)

    query_profiles = {}
    for ct in query_celltypes:
        query_mask = (query_labels == ct)
        if np.sum(query_mask) > 0:
            if sparse.issparse(X_query):
                query_profiles[ct] = np.array(X_query[query_mask].mean(axis=0)).flatten()
            else:
                query_profiles[ct] = X_query[query_mask].mean(axis=0)

    for i, ref_ct in enumerate(ref_celltypes):
        if ref_ct not in ref_profiles:
            continue

        for j, query_ct in enumerate(query_celltypes):
            if query_ct not in query_profiles:
                continue

            ref_vec = ref_profiles[ref_ct].reshape(1, -1)
            query_vec = query_profiles[query_ct].reshape(1, -1)

            similarity = cosine_similarity(ref_vec, query_vec)[0, 0]
            similarity_matrix[i, j] = similarity

    return ref_celltypes, query_celltypes, similarity_matrix

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--sample1", required=True, help="Reference sample name")
    parser.add_argument("--sample2", required=True, help="Query sample name")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    celltype_col = "celltype"
    sample_col = "sample"

    adata_ref = adata[adata.obs[sample_col] == args.sample1].copy()
    adata_query = adata[adata.obs[sample_col] == args.sample2].copy()

    if adata_ref.n_obs == 0 or adata_query.n_obs == 0:
        print("ERROR: empty reference or query sample")
        sys.exit(1)

    print(f"Reference ({args.sample1}): {adata_ref.shape}")
    print(f"Query ({args.sample2}): {adata_query.shape}")

    print("Identifying highly variable genes...")
    
    adata_combined = sc.concat([adata_ref, adata_query], join='inner', label='sample_origin')
    
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor='seurat')
    
    hvg_mask = adata_combined.var['highly_variable']
    hvg_names = adata_combined.var_names[hvg_mask].tolist()
    
    print(f"Number of highly variable genes: {len(hvg_names)}")
    
    adata_ref = adata_ref[:, hvg_names].copy()
    adata_query = adata_query[:, hvg_names].copy()
    
    print(f"After HVG selection - Ref: {adata_ref.shape}, Query: {adata_query.shape}")

    X_ref = adata_ref.X
    X_query = adata_query.X

    ref_labels = adata_ref.obs[celltype_col].values
    query_labels = adata_query.obs[celltype_col].values

    print("Calculating cosine similarity matrix...")
    ref_cts, query_cts, sim_matrix = calculate_cosine_similarity_sparse(X_ref, X_query, ref_labels, query_labels)

    sim_df = pd.DataFrame(sim_matrix, index=ref_cts, columns=query_cts)

    print("Plotting similarity matrix...")
    plt.figure(figsize=(max(12, len(sim_df.columns) * 0.8), max(10, len(sim_df.index) * 0.7)))

    plot_data = sim_df.fillna(0).values

    im = plt.imshow(plot_data, aspect='auto', cmap='viridis')
    plt.colorbar(im, label="Cosine Similarity Score", shrink=0.8)

    plt.xticks(range(len(sim_df.columns)), sim_df.columns, rotation=45, ha='right', fontsize=10)
    plt.yticks(range(len(sim_df.index)), sim_df.index, fontsize=10)

    plt.xlabel(f"Query: {args.sample2}", fontsize=12, fontweight='bold')
    plt.ylabel(f"Reference: {args.sample1}", fontsize=12, fontweight='bold')
    plt.title(f"Cell Type Similarity: {args.sample1} → {args.sample2}\n(Cosine Similarity - Highly Variable Genes)", fontsize=14, fontweight='bold', pad=20)

    for i in range(len(sim_df.index)):
        for j in range(len(sim_df.columns)):
            val = plot_data[i, j]
            plt.text(j, i, f"{val:.3f}",
                     ha="center", va="center",
                     color="black",
                     fontsize=8 if len(sim_df.index) < 20 else 6,
                     fontweight='bold' if val > 0.7 else 'normal')

    plt.tight_layout()
    
    os.makedirs("figures", exist_ok=True)
    plt.savefig(f"figures/{args.prefix}_cosine_similarity.png", dpi=300, bbox_inches='tight')
    plt.close()

    csv_path = f"{args.prefix}_cosine_similarity.csv"
    sim_df.to_csv(csv_path)
    print(f"Saved similarity matrix to: {csv_path}")
    print(f"Saved heatmap to: figures/{args.prefix}_cosine_similarity.png")

if __name__ == "__main__":
    main()
