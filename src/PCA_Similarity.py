import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)

LAYER = "log1p"

def get_X(adata):
    X = adata.layers[LAYER] if LAYER in adata.layers else adata.X
    return X.toarray() if issparse(X) else X

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample1", required=True, help="Name of first sample")
    parser.add_argument("--sample2", required=True, help="Name of second sample")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCA components")
    parser.add_argument("--max_cells", type=int, default=10000, help="Max cells per cell type")
    args = parser.parse_args()

    print("Reading input file...", flush=True)
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes", flush=True)

    if "celltype" not in adata.obs:
        raise ValueError("Missing celltype column")
    if "sample" not in adata.obs:
        raise ValueError("Missing sample column")

    # Extract the two samples
    print(f"Filtering samples...", flush=True)
    adata1 = adata[adata.obs["sample"] == args.sample1].copy()
    adata2 = adata[adata.obs["sample"] == args.sample2].copy()
    print(f"  {args.sample1}: {adata1.n_obs} cells", flush=True)
    print(f"  {args.sample2}: {adata2.n_obs} cells", flush=True)

    # Get common cell types
    celltypes1 = set(adata1.obs["celltype"].unique())
    celltypes2 = set(adata2.obs["celltype"].unique())
    common_celltypes = sorted(list(celltypes1.intersection(celltypes2)))
    print(f"Common cell types: {common_celltypes}", flush=True)

    # ------------------------------------------------
    # 1. Build combined matrix for PCA (no concat)
    # ------------------------------------------------
    print("\nBuilding combined matrix for PCA...", flush=True)
    rng = np.random.default_rng(42)
    
    all_X = []
    for ct in common_celltypes:
        sub1 = adata1[adata1.obs["celltype"] == ct]
        sub2 = adata2[adata2.obs["celltype"] == ct]
        
        # Subsample each group
        n1 = min(sub1.n_obs, args.max_cells)
        n2 = min(sub2.n_obs, args.max_cells)
        if n1 < sub1.n_obs:
            sub1 = sub1[rng.choice(sub1.n_obs, n1, replace=False)]
        if n2 < sub2.n_obs:
            sub2 = sub2[rng.choice(sub2.n_obs, n2, replace=False)]
        
        X1 = get_X(sub1)
        X2 = get_X(sub2)
        all_X.append(X1)
        all_X.append(X2)
    
    X_combined = np.vstack(all_X)
    print(f"  Combined matrix shape: {X_combined.shape[0]} cells, {X_combined.shape[1]} genes", flush=True)
    
    # Remove near-constant genes
    var = np.var(X_combined, axis=0)
    keep = var > 1e-8
    X_combined = X_combined[:, keep]
    print(f"  Genes after variance filter: {X_combined.shape[1]}", flush=True)
    
    # Standardise and fit PCA once
    print("Fitting PCA...", flush=True)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_combined)
    pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    X_pca = pca.fit_transform(X_scaled)
    print(f"  PCA fitted with {pca.n_components_} components", flush=True)

    # ------------------------------------------------
    # 2. Project all cells into this fixed PCA space
    # ------------------------------------------------
    def project_group(adata_group):
        X = get_X(adata_group)
        X = X[:, keep]
        X_scaled = scaler.transform(X)
        return pca.transform(X_scaled)

    # Store projections
    proj1 = {}
    proj2 = {}
    
    print("\nProjecting Sample1 cells...", flush=True)
    for ct in common_celltypes:
        sub = adata1[adata1.obs["celltype"] == ct]
        if sub.n_obs > args.max_cells:
            sub = sub[rng.choice(sub.n_obs, args.max_cells, replace=False)]
        proj1[ct] = project_group(sub)
        print(f"  {ct}: {proj1[ct].shape[0]} cells, {proj1[ct].shape[1]} PCs", flush=True)

    print("\nProjecting Sample2 cells...", flush=True)
    for ct in common_celltypes:
        sub = adata2[adata2.obs["celltype"] == ct]
        if sub.n_obs > args.max_cells:
            sub = sub[rng.choice(sub.n_obs, args.max_cells, replace=False)]
        proj2[ct] = project_group(sub)
        print(f"  {ct}: {proj2[ct].shape[0]} cells, {proj2[ct].shape[1]} PCs", flush=True)

    # ------------------------------------------------
    # 3. Calculate pairwise distances
    # ------------------------------------------------
    def calc_wdist(projA, projB):
        n_pcs = min(projA.shape[1], projB.shape[1])
        dists = [wasserstein_distance(projA[:, i], projB[:, i]) for i in range(n_pcs)]
        return np.mean(dists)

    n_ct = len(common_celltypes)
    dist_matrix = np.zeros((n_ct, n_ct))
    
    print("\nCalculating pairwise distances...", flush=True)
    for i, ct1 in enumerate(common_celltypes):
        for j, ct2 in enumerate(common_celltypes):
            dist_matrix[i, j] = calc_wdist(proj1[ct1], proj2[ct2])
        print(f"  Done: {ct1}", flush=True)

    # FIX 1: Make the distance matrix symmetric
    dist_matrix = (dist_matrix + dist_matrix.T) / 2

    # Normalise to similarity
    max_dist = dist_matrix.max()
    if max_dist > 0:
        sim_matrix = 1 - (dist_matrix / max_dist)
    else:
        sim_matrix = np.ones_like(dist_matrix)

    # ------------------------------------------------
    # 4. Plots
    # ------------------------------------------------
    os.makedirs("figures", exist_ok=True)
    
    print("\nGenerating plots...", flush=True)
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(max(10, n_ct * 0.6), max(8, n_ct * 0.5)))
    im = ax.imshow(sim_matrix, cmap='viridis', aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(n_ct))
    ax.set_yticks(range(n_ct))
    ax.set_xticklabels(common_celltypes, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(common_celltypes, fontsize=9)
    for i in range(n_ct):
        for j in range(n_ct):
            text_color = "white" if sim_matrix[i, j] < 0.5 else "black"
            ax.text(j, i, f'{sim_matrix[i, j]:.2f}', ha="center", va="center", color=text_color, fontsize=8)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Similarity Score (1 = identical, 0 = completely different)', fontsize=10)
    plt.title(f'Cell Type Transcriptional Similarity\n{args.prefix}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    sim_heatmap = os.path.join("figures", f"{args.prefix}_celltype_similarity_heatmap.png")
    plt.savefig(sim_heatmap, dpi=300, bbox_inches='tight')
    print(f"  Saved: {sim_heatmap}", flush=True)
    plt.close()

    # PCA visualization (first 2 PCs)
    all_pcs = []
    all_labels = []
    for ct in common_celltypes:
        p1 = proj1[ct][:, :2]
        p2 = proj2[ct][:, :2]
        all_pcs.append(p1)
        all_pcs.append(p2)
        all_labels.extend([f"{ct}_{args.sample1}"] * len(p1))
        all_labels.extend([f"{ct}_{args.sample2}"] * len(p2))
    all_pcs = np.vstack(all_pcs)

    fig, ax = plt.subplots(figsize=(10, 8))
    unique_labels = list(set(all_labels))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))
    for i, label in enumerate(unique_labels):
        mask = np.array(all_labels) == label
        ax.scatter(all_pcs[mask, 0], all_pcs[mask, 1], c=[colors[i]], label=label, alpha=0.6, s=10)
    ax.set_xlabel('PC1', fontsize=10)
    ax.set_ylabel('PC2', fontsize=10)
    ax.set_title(f'PCA Projection - All Cell Types\n{args.prefix}', fontsize=12, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    pca_file = os.path.join("figures", f"{args.prefix}_pca_all_celltypes.png")
    plt.savefig(pca_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {pca_file}", flush=True)
    plt.close()

    # FIX 2: Diagonal bar plot with wider figure to prevent x-label overlap
    fig, ax = plt.subplots(figsize=(max(12, n_ct * 0.8), 6))
    diag_sims = [sim_matrix[i, i] for i in range(n_ct)]
    bars = ax.bar(common_celltypes, diag_sims, color=plt.cm.viridis(np.linspace(0.2, 0.8, n_ct)))
    ax.set_xlabel('Cell Type', fontsize=10)
    ax.set_ylabel('Cross-Sample Similarity', fontsize=10)
    ax.set_title(f'Same Cell Type Across Samples\n{args.sample1} vs {args.sample2}', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1.05])
    ax.tick_params(axis='x', rotation=45, labelsize=10)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
    ax.grid(alpha=0.3, axis='y')
    for bar, val in zip(bars, diag_sims):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.3f}', ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    bar_file = os.path.join("figures", f"{args.prefix}_diagonal_similarities.png")
    plt.savefig(bar_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {bar_file}", flush=True)
    plt.close()

    # Summary
    print("\n" + "="*50, flush=True)
    print(f"Same Cell Type Similarity: {args.sample1} vs {args.sample2}", flush=True)
    print("="*50, flush=True)
    for i, ct in enumerate(common_celltypes):
        print(f"{ct:20s}: {sim_matrix[i, i]:.4f}", flush=True)
    print("="*50, flush=True)
    print(f"\nAverage: {np.mean(diag_sims):.4f} ± {np.std(diag_sims):.4f}", flush=True)
    
    print("\nDone!", flush=True)

if __name__ == "__main__":
    main()
