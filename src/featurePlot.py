#!/usr/bin/env python3

import scanpy as sc
import matplotlib.pyplot as plt
import argparse
import os


def load_markers(path):
    with open(path, "r") as f:
        return [g.strip() for g in f if g.strip()]


def plot_gene_two_sample_panels(adata, gene, outdir, prefix):
    if gene not in adata.var_names:
        print(f"[WARN] {gene} not in dataset")
        return

    if "sample" not in adata.obs:
        raise ValueError("adata.obs['sample'] is missing")

    samples = list(adata.obs["sample"].unique())

    if len(samples) != 2:
        print(f"[WARN] Expected 2 samples, got {len(samples)}")

    fig, axes = plt.subplots(1, len(samples), figsize=(12, 5))

    if len(samples) == 1:
        axes = [axes]

    # FIX: global scaling across BOTH samples
    gene_data = adata[:, gene].X
    vmin = gene_data.min()
    vmax = gene_data.max()

    for i, s in enumerate(samples):
        ax = axes[i]

        ad = adata[adata.obs["sample"] == s]

        sc.pl.umap(
            ad,
            color=gene,
            ax=ax,
            show=False,
            frameon=False,
            title=str(s),
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            size=20
        )

    fig.suptitle(f"{prefix} | {gene}")
    plt.tight_layout()

    outpath = os.path.join(outdir, f"{prefix}_{gene}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--markers", required=True)

    args = parser.parse_args()

    print("[INFO] loading h5ad")
    adata = sc.read_h5ad(args.input)

    print("[INFO] loading markers")
    genes = load_markers(args.markers)

    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)

    if "X_umap" not in adata.obsm:
        print("[INFO] computing UMAP")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    for g in genes:
        print(f"[INFO] plotting {g}")
        plot_gene_two_sample_panels(adata, g, outdir, args.prefix)

    print("[DONE]")


if __name__ == "__main__":
    main()
