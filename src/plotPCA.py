#!/usr/bin/env python3
import scanpy as sc
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    outdir = os.path.abspath("figures")
    os.makedirs(outdir, exist_ok=True)

    if "PCs" not in adata.varm:
        raise ValueError("Missing adata.varm['PCs']")

    loadings = adata.varm["PCs"]

    if loadings.shape[1] < 3:
        raise ValueError("Need at least 3 PCs")

    df = pd.DataFrame({
        "PC1": loadings[:, 0],
        "PC2": loadings[:, 1],
        "PC3": loadings[:, 2],
        "gene": adata.var_names
    })

    df["abs_PC1"] = np.abs(df["PC1"])
    df["abs_PC2"] = np.abs(df["PC2"])
    df["abs_PC3"] = np.abs(df["PC3"])

    pc1_top = df.nlargest(10, "abs_PC1")
    pc2_top = df.nlargest(10, "abs_PC2")
    pc3_top = df.nlargest(10, "abs_PC3")

    pc1_pos = df.nlargest(10, "PC1")
    pc1_neg = df.nsmallest(10, "PC1")

    pc2_pos = df.nlargest(10, "PC2")
    pc2_neg = df.nsmallest(10, "PC2")

    pc3_pos = df.nlargest(10, "PC3")
    pc3_neg = df.nsmallest(10, "PC3")

    print("\nPC1 top drivers (absolute):")
    print(list(pc1_top["gene"]))
    print("PC1 positive genes:")
    print(list(pc1_pos["gene"]))
    print("PC1 negative genes:")
    print(list(pc1_neg["gene"]))

    print("\nPC2 top drivers (absolute):")
    print(list(pc2_top["gene"]))
    print("PC2 positive genes:")
    print(list(pc2_pos["gene"]))
    print("PC2 negative genes:")
    print(list(pc2_neg["gene"]))

    print("\nPC3 top drivers (absolute):")
    print(list(pc3_top["gene"]))
    print("PC3 positive genes:")
    print(list(pc3_pos["gene"]))
    print("PC3 negative genes:")
    print(list(pc3_neg["gene"]))

    def save(fig, name):
        path = os.path.join(outdir, f"{args.prefix}_{name}.png")
        fig.tight_layout()
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("Saved:", path)

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].barh(pc1_pos["gene"], pc1_pos["PC1"])
    ax[0].invert_yaxis()
    ax[1].barh(pc1_neg["gene"], pc1_neg["PC1"])
    ax[1].invert_yaxis()
    save(fig, "pc1_genes")

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].barh(pc2_pos["gene"], pc2_pos["PC2"])
    ax[0].invert_yaxis()
    ax[1].barh(pc2_neg["gene"], pc2_neg["PC2"])
    ax[1].invert_yaxis()
    save(fig, "pc2_genes")

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].barh(pc3_pos["gene"], pc3_pos["PC3"])
    ax[0].invert_yaxis()
    ax[1].barh(pc3_neg["gene"], pc3_neg["PC3"])
    ax[1].invert_yaxis()
    save(fig, "pc3_genes")

    if "batch" in adata.obs:
        sc.pl.pca(
            adata,
            color="batch",
            components=["1,2"],
            show=False,
            legend_loc="right margin"
        )
        plt.savefig(os.path.join(outdir, f"{args.prefix}_pca12_batch.png"),
                    dpi=120, bbox_inches="tight")
        plt.close()

        sc.pl.pca(
            adata,
            color="batch",
            components=["1,3"],
            show=False,
            legend_loc="right margin"
        )
        plt.savefig(os.path.join(outdir, f"{args.prefix}_pca13_batch.png"),
                    dpi=120, bbox_inches="tight")
        plt.close()

        # ✅ ONLY ADDITION (PC1 vs PC2 batch split view)
        pc1 = adata.obsm["X_pca"][:, 0]
        pc2 = adata.obsm["X_pca"][:, 1]

        plt.figure(figsize=(6, 5))
        for b in adata.obs["batch"].unique():
            idx = adata.obs["batch"] == b
            plt.scatter(pc1[idx], pc2[idx], label=b, s=5)

        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.legend()
        plt.title("PC1 vs PC2 colored by batch")

        plt.savefig(os.path.join(outdir, f"{args.prefix}_pc1_pc2_batch.png"),
                    dpi=120, bbox_inches="tight")
        plt.close()

        print("Saved PC1 vs PC2 batch split plot")

    print("DONE — ALL FILES SAVED IN:", outdir)


if __name__ == "__main__":
    main()
