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


def compute_similarity_matrix(mat, method="cosine"):
    if method == "cosine":
        return cosine_similarity(mat, mat)

    elif method == "pearson":
        return np.corrcoef(mat)

    elif method == "spearman":
        ranked = np.apply_along_axis(
            lambda x: pd.Series(x).rank().values,
            1,
            mat
        )
        return np.corrcoef(ranked)


def compute_state_pseudobulk(adata, celltype_col, sample_col):
    states = adata.obs[celltype_col].astype(str) + "_" + adata.obs[sample_col].astype(str)
    adata.obs["state"] = states

    profiles = {}

    for state in adata.obs["state"].unique():
        mask = adata.obs["state"] == state
        X = adata[mask].X

        if sparse.issparse(X):
            profiles[state] = np.array(X.mean(axis=0)).flatten()
        else:
            profiles[state] = X.mean(axis=0)

    return profiles


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--method", default="cosine",
                        choices=["cosine", "pearson", "spearman"])
    parser.add_argument("--hvg", type=int, default=0, 
                        help="Number of highly variable genes to use (0 = use all genes)")

    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    celltype_col = "celltype"
    sample_col = "sample"

    # ----------------------------
    # OPTIONAL: FILTER TO HIGHLY VARIABLE GENES
    # ----------------------------
    if args.hvg > 0:
        print(f"Selecting top {args.hvg} HVGs...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.hvg,
            flavor="seurat"
        )
        genes = adata.var_names[adata.var["highly_variable"]]
        adata = adata[:, genes].copy()
        print(f"Kept {adata.shape[1]} genes")
    else:
        print("Using all genes")

    print("Computing state pseudobulk...")
    profiles = compute_state_pseudobulk(adata, celltype_col, sample_col)

    states = list(profiles.keys())
    mat = np.vstack([profiles[s] for s in states])

    print(f"Computing {args.method} similarity...")

    sim_matrix = compute_similarity_matrix(mat, args.method)

    sim_df = pd.DataFrame(
        sim_matrix,
        index=states,
        columns=states
    )

    # -------------------------
    # SAVE CSV
    # -------------------------
    hvg_text = f"_hvg{args.hvg}" if args.hvg > 0 else ""
    csv_path = f"{args.prefix}_{args.method}{hvg_text}_similarity.csv"
    sim_df.to_csv(csv_path)

    print(f"Saved: {csv_path}")

    # -------------------------
    # PLOT
    # -------------------------
    plt.figure(figsize=(
        max(10, len(states) * 0.6),
        max(10, len(states) * 0.6)
    ))

    data = sim_df.values

    im = plt.imshow(data, cmap="viridis", aspect="auto")
    plt.colorbar(label=f"{args.method} similarity")

    plt.xticks(range(len(states)), states, rotation=90, fontsize=6)
    plt.yticks(range(len(states)), states, fontsize=6)

    # -------------------------
    # ADD NUMBERS IN CELLS
    # -------------------------
    for i in range(len(states)):
        for j in range(len(states)):
            plt.text(
                j, i,
                f"{data[i, j]:.2f}",
                ha="center",
                va="center",
                color="black",
                fontsize=5
            )

    hvg_title = f"HVG={args.hvg}" if args.hvg > 0 else "All genes"
    plt.title(f"Celltype × Sample State Similarity ({args.method}) - {hvg_title}")

    plt.tight_layout()

    os.makedirs("figures", exist_ok=True)

    plt.savefig(
        f"figures/{args.prefix}_{args.method}{hvg_text}_similarity.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Done.")


if __name__ == "__main__":
    main()
