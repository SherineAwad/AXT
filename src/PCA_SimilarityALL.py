#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import pairwise_kernels
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)


# ----------------------------
# GET MATRIX
# ----------------------------
def get_X(adata):
    X = adata.X
    return X.toarray() if issparse(X) else X


# ----------------------------
# MMD (simple kernel version)
# ----------------------------
def compute_mmd(X, Y):
    # RBF kernel (default gamma auto)
    Kxx = pairwise_kernels(X, X, metric="rbf")
    Kyy = pairwise_kernels(Y, Y, metric="rbf")
    Kxy = pairwise_kernels(X, Y, metric="rbf")

    return Kxx.mean() + Kyy.mean() - 2 * Kxy.mean()


# ----------------------------
# MAIN
# ----------------------------
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--method", choices=["wasserstein", "mmd"], default="wasserstein")
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--max_cells", type=int, default=10000)
    parser.add_argument("--hvg", type=int, default=0, help="Number of highly variable genes to use (0 = use all genes)")

    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    assert "celltype" in adata.obs
    assert "sample" in adata.obs

    # ----------------------------
    # OPTIONAL: FILTER TO HIGHLY VARIABLE GENES
    # ----------------------------
    if args.hvg > 0:
        print(f"Filtering to top {args.hvg} highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, flavor='seurat')
        adata = adata[:, adata.var.highly_variable]
        print(f"Kept {adata.shape[1]} genes")
    else:
        print("Using all genes")

    # ----------------------------
    # BUILD STATES
    # ----------------------------
    adata.obs["state"] = (
        adata.obs["celltype"].astype(str)
        + "_"
        + adata.obs["sample"].astype(str)
    )

    states = sorted(adata.obs["state"].unique())
    print(f"Found {len(states)} states")

    rng = np.random.default_rng(42)

    # ----------------------------
    # COLLECT DATA
    # ----------------------------
    state_data = {}
    all_X = []

    for st in states:
        sub = adata[adata.obs["state"] == st]

        if sub.n_obs > args.max_cells:
            idx = rng.choice(sub.n_obs, args.max_cells, replace=False)
            sub = sub[idx]

        X = get_X(sub)
        state_data[st] = X
        all_X.append(X)

    X_all = np.vstack(all_X)

    print("Fitting PCA...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_all)

    pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    pca.fit(X_scaled)

    # ----------------------------
    # PROJECT STATES
    # ----------------------------
    proj = {}

    for st in states:
        X = state_data[st]
        Xs = scaler.transform(X)
        proj[st] = pca.transform(Xs)

    # ----------------------------
    # DISTANCE MATRIX
    # ----------------------------
    n = len(states)
    dist = np.zeros((n, n))

    print(f"Computing {args.method} distances...")

    for i, s1 in enumerate(states):
        for j, s2 in enumerate(states):

            A = proj[s1]
            B = proj[s2]

            if args.method == "wasserstein":
                n_pcs = min(A.shape[1], B.shape[1])
                dists = [
                    wasserstein_distance(A[:, k], B[:, k])
                    for k in range(n_pcs)
                ]
                dist[i, j] = np.mean(dists)

            elif args.method == "mmd":
                dist[i, j] = compute_mmd(A, B)

        print(f"done: {s1}")

    # ----------------------------
    # SYMMETRIZE
    # ----------------------------
    dist = (dist + dist.T) / 2

    # ----------------------------
    # NORMALISE TO SIMILARITY
    # ----------------------------
    max_d = dist.max()
    sim = 1 - dist / max_d if max_d > 0 else np.ones_like(dist)

    # ----------------------------
    # PLOT
    # ----------------------------
    os.makedirs("figures", exist_ok=True)

    plt.figure(figsize=(max(10, n * 0.6), max(10, n * 0.6)))

    plt.imshow(sim, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    plt.colorbar(label=f"{args.method} similarity")

    plt.xticks(range(n), states, rotation=90, fontsize=6)
    plt.yticks(range(n), states, fontsize=6)

    # ----------------------------
    # CELL VALUES (rounded)
    # ----------------------------
    for i in range(n):
        for j in range(n):
            plt.text(
                j, i,
                f"{sim[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=5,
                color="black"
            )

    hvg_text = f"_hvg{args.hvg}" if args.hvg > 0 else ""
    plt.title(f"PCA + {args.method.upper()} State Similarity (HVG={args.hvg if args.hvg > 0 else 'all'})\n{args.prefix}")

    plt.tight_layout()

    out = f"figures/{args.prefix}_pca_{args.method}{hvg_text}.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
