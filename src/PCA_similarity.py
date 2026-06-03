#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_kernels, pairwise_distances
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)


# ----------------------------
# GET MATRIX
# ----------------------------
def get_X(adata):
    X = adata.X
    if issparse(X):
        return X.toarray()
    return np.asarray(X)


# ----------------------------
# MMD ONLY
# ----------------------------
def compute_mmd(X, Y):
    Z = np.vstack([X, Y])
    dists = pairwise_distances(Z, metric="euclidean")
    gamma = 1.0 / (np.median(dists) ** 2 + 1e-8)

    Kxx = pairwise_kernels(X, X, metric="rbf", gamma=gamma)
    Kyy = pairwise_kernels(Y, Y, metric="rbf", gamma=gamma)
    Kxy = pairwise_kernels(X, Y, metric="rbf", gamma=gamma)

    return Kxx.mean() + Kyy.mean() - 2 * Kxy.mean()


# ----------------------------
# MAIN
# ----------------------------
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--max_cells", type=int, default=10000)
    parser.add_argument("--hvg", type=int, default=0)

    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    assert "celltype" in adata.obs
    assert "sample" in adata.obs

    if args.hvg > 0:
        sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, flavor="seurat")
        adata = adata[:, adata.var.highly_variable]

    adata.obs["state"] = (
        adata.obs["celltype"].astype(str)
        + "_"
        + adata.obs["sample"].astype(str)
    )

    states = sorted(adata.obs["state"].unique())
    rng = np.random.default_rng(42)

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

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_all)

    pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    pca.fit(X_scaled)

    proj = {}

    for st in states:
        X = state_data[st]
        Xs = scaler.transform(X)
        proj[st] = pca.transform(Xs)

    n = len(states)
    dist = np.zeros((n, n))

    for i, s1 in enumerate(states):
        for j in range(i, n):

            s2 = states[j]

            A = proj[s1]
            B = proj[s2]

            d = compute_mmd(A, B)

            dist[i, j] = d
            dist[j, i] = d

        print(f"done: {s1}")

    sim = np.exp(-dist)

    os.makedirs("figures", exist_ok=True)

    plt.figure(figsize=(max(10, n * 0.6), max(10, n * 0.6)))

    plt.imshow(sim, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    plt.colorbar(label="MMD similarity")

    plt.xticks(range(n), states, rotation=90, fontsize=6)
    plt.yticks(range(n), states, fontsize=6)

    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{sim[i, j]:.2f}",
                     ha="center", va="center", fontsize=5, color="black")

    plt.title(f"PCA + MMD State Similarity\n{args.prefix}")

    plt.tight_layout()

    out = f"figures/{args.prefix}_pca_mmd.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
