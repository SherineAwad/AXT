#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import warnings
import ot

warnings.filterwarnings("ignore")


def get_X(adata, layer):
    X = adata.layers[layer] if layer in adata.layers else adata.X
    return X.toarray() if issparse(X) else X


def pot_score(A, B, max_cells, min_cells):
    try:
        if A is None or B is None:
            return np.nan
        if len(A) < min_cells or len(B) < min_cells:
            return np.nan
        
        if len(A) > max_cells:
            idx = np.random.choice(len(A), max_cells, replace=False)
            A = A[idx]
        if len(B) > max_cells:
            idx = np.random.choice(len(B), max_cells, replace=False)
            B = B[idx]
        
        scaler = StandardScaler()
        A_norm = scaler.fit_transform(A)
        B_norm = scaler.transform(B)
        
        M = ot.dist(A_norm, B_norm, metric='euclidean')
        M = M / M.max()
        
        a = np.ones(len(A_norm)) / len(A_norm)
        b = np.ones(len(B_norm)) / len(B_norm)
        
        T = ot.sinkhorn(a, b, M, reg=0.1)
        
        cost = np.sum(T * M)
        similarity = 1 - cost
        
        return max(0.0, min(1.0, similarity))
        
    except Exception:
        return np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument("--layer", default="log1p", help="Layer to use for expression")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCs")
    parser.add_argument("--max_cells", type=int, default=10000, help="Maximum cells to subsample per group")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum cells required per group")
    parser.add_argument("--hvg", type=int, default=0, help="Number of highly variable genes to use (0 = use all genes)")
    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    # =====================================================
    # OPTIONAL: FILTER TO HIGHLY VARIABLE GENES
    # =====================================================
    if args.hvg > 0:
        print(f"Filtering to top {args.hvg} highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, batch_key='sample', flavor='seurat')
        adata = adata[:, adata.var.highly_variable]
        print(f"Kept {adata.shape[1]} genes")
    else:
        print("Using all genes")

    # =====================================================
    # FILTER OUT ZERO-CELL COMBINATIONS
    # =====================================================
    print("Filtering out celltype-sample combinations with zero cells...")
    
    valid_groups = adata.obs.groupby(['celltype', 'sample']).size()
    valid_groups = valid_groups[valid_groups > 0]
    
    print("Keeping these combinations:")
    for (ct, s), n in valid_groups.items():
        print(f"  {ct}.{s}: {n} cells")
    
    celltypes = sorted(set([ct for ct, s in valid_groups.index]))
    samples = sorted(set([s for ct, s in valid_groups.index]))
    
    print(f"\nCelltypes with cells: {celltypes}")
    print(f"Samples with cells: {samples}")

    # =====================================================
    # BUILD PCA SPACE
    # =====================================================
    print("\nBuilding global PCA space...")
    
    all_cells = []
    for ct in celltypes:
        for s in samples:
            sub = adata[(adata.obs["celltype"] == ct) & (adata.obs["sample"] == s)]
            if sub.n_obs > 0:
                all_cells.append(get_X(sub, args.layer))
    
    X_all = np.vstack(all_cells)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_all)
    pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    X_pca = pca.fit_transform(X_scaled)
    print(f"PCA shape: {X_pca.shape}")

    # =====================================================
    # CREATE STATES DICTIONARY
    # =====================================================
    states = {}
    idx = 0
    for ct in celltypes:
        for s in samples:
            sub = adata[(adata.obs["celltype"] == ct) & (adata.obs["sample"] == s)]
            n = sub.n_obs
            if n > 0:
                states[(ct, s)] = X_pca[idx:idx + n]
                idx += n
    
    keys = list(states.keys())
    n_states = len(keys)
    print(f"\nTotal states to compare: {n_states}")
    
    # =====================================================
    # COMPUTE SIMILARITY MATRIX
    # =====================================================
    matrix = np.zeros((n_states, n_states))
    
    print("\nComputing POT similarity...")
    print(f"Parameters: max_cells={args.max_cells}, min_cells={args.min_cells}, hvg={args.hvg}\n")
    
    for i, k1 in enumerate(keys):
        for j, k2 in enumerate(keys):
            if i == j:
                matrix[i, j] = 1.0
            elif i < j:
                matrix[i, j] = pot_score(states[k1], states[k2], args.max_cells, args.min_cells)
                matrix[j, i] = matrix[i, j]
        print(f"  done: {k1[0]}.{k1[1]}")
    
    # =====================================================
    # PLOT
    # =====================================================
    os.makedirs("figures", exist_ok=True)
    
    labels = [f"{ct}.{s}" for (ct, s) in keys]
    n = len(labels)
    
    plt.figure(figsize=(max(8, n * 0.4), max(8, n * 0.4)))
    plt.imshow(matrix, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="POT similarity")
    plt.xticks(range(n), labels, rotation=90, fontsize=8)
    plt.yticks(range(n), labels, fontsize=8)
    
    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            if not np.isnan(val):
                color = "white" if val < 0.5 else "black"
                plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color=color)
    
    hvg_text = f"_hvg{args.hvg}" if args.hvg > 0 else ""
    plt.title(f"POT similarity (HVG={args.hvg if args.hvg > 0 else 'all'})\n{args.prefix}")
    plt.tight_layout()
    
    out = f"figures/{args.prefix}_pot{hvg_text}_matrix.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"\nSaved: {out}")
    print(f"Matrix shape: {n_states} x {n_states}")
    print("Done.")


if __name__ == "__main__":
    main()
