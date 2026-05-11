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

def get_pca_projection(adata, n_pcs, max_cells):
    """Helper function to compute PCA projection for a cell type"""
    if adata.n_obs > max_cells:
        rng = np.random.default_rng(42)
        sample_idx = rng.choice(adata.n_obs, max_cells, replace=False)
        adata = adata[sample_idx].copy()
    
    X = get_X(adata)
    var = np.var(X, axis=0)
    keep = var > 1e-8
    X = X[:, keep]
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=min(n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    X_pca = pca.fit_transform(X_scaled)
    
    return X_pca

def calculate_wasserstein_distance(proj1, proj2):
    n_pcs = min(proj1.shape[1], proj2.shape[1])
    distances = [wasserstein_distance(proj1[:, pc], proj2[:, pc]) for pc in range(n_pcs)]
    return np.mean(distances)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample1", required=True, help="Name of first sample")
    parser.add_argument("--sample2", required=True, help="Name of second sample")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCA components")
    parser.add_argument("--max_cells", type=int, default=10000, help="Max cells per cell type")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if "celltype" not in adata.obs:
        raise ValueError("Missing celltype column")
    if "sample" not in adata.obs:
        raise ValueError("Missing sample column")

    # Extract the two samples
    adata1 = adata[adata.obs["sample"] == args.sample1].copy()
    adata2 = adata[adata.obs["sample"] == args.sample2].copy()

    print(f"Sample1 ({args.sample1}): {adata1.n_obs} cells")
    print(f"Sample2 ({args.sample2}): {adata2.n_obs} cells")

    # Get common cell types
    celltypes1 = set(adata1.obs["celltype"].unique())
    celltypes2 = set(adata2.obs["celltype"].unique())
    common_celltypes = sorted(list(celltypes1.intersection(celltypes2)))
    
    print(f"Common cell types: {common_celltypes}")

    # Store PCA projections
    pca_projections_sample1 = {}
    pca_projections_sample2 = {}

    print("\nProcessing Sample1:")
    for ct in common_celltypes:
        print(f"  {ct}")
        sub = adata1[adata1.obs["celltype"] == ct].copy()
        pca_projections_sample1[ct] = get_pca_projection(sub, args.n_pcs, args.max_cells)
        print(f"    {pca_projections_sample1[ct].shape[0]} cells, {pca_projections_sample1[ct].shape[1]} PCs kept")

    print("\nProcessing Sample2:")
    for ct in common_celltypes:
        print(f"  {ct}")
        sub = adata2[adata2.obs["celltype"] == ct].copy()
        pca_projections_sample2[ct] = get_pca_projection(sub, args.n_pcs, args.max_cells)
        print(f"    {pca_projections_sample2[ct].shape[0]} cells, {pca_projections_sample2[ct].shape[1]} PCs kept")

    # Calculate pairwise distances
    n_ct = len(common_celltypes)
    distance_matrix = np.zeros((n_ct, n_ct))
    
    print("\nCalculating pairwise distances...")
    for i, ct1 in enumerate(common_celltypes):
        for j, ct2 in enumerate(common_celltypes):
            proj1 = pca_projections_sample1[ct1]
            proj2 = pca_projections_sample2[ct2]
            distance_matrix[i, j] = calculate_wasserstein_distance(proj1, proj2)
        print(f"  Done: {ct1}")

    # Normalize distances
    max_dist = distance_matrix.max()
    if max_dist > 0:
        similarity_matrix = 1 - (distance_matrix / max_dist)
    else:
        similarity_matrix = np.ones_like(distance_matrix)

    # -------------------------
    # PLOTS
    # -------------------------
    os.makedirs("figures", exist_ok=True)
    
    # 1. Heatmap
    fig, ax = plt.subplots(figsize=(max(10, n_ct * 0.6), max(8, n_ct * 0.5)))
    im = ax.imshow(similarity_matrix, cmap='viridis', aspect='auto', vmin=0, vmax=1)
    
    ax.set_xticks(range(n_ct))
    ax.set_yticks(range(n_ct))
    ax.set_xticklabels(common_celltypes, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(common_celltypes, fontsize=9)
    
    for i in range(n_ct):
        for j in range(n_ct):
            text_color = "white" if similarity_matrix[i, j] < 0.5 else "black"
            ax.text(j, i, f'{similarity_matrix[i, j]:.2f}',
                   ha="center", va="center", color=text_color, fontsize=8)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Similarity Score (1 = identical, 0 = completely different)', fontsize=10)
    plt.title(f'Cell Type Transcriptional Similarity\n{args.prefix}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    sim_heatmap = os.path.join("figures", f"{args.prefix}_celltype_similarity_heatmap.png")
    plt.savefig(sim_heatmap, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {sim_heatmap}")
    plt.close()
    
    # 2. PCA visualization
    all_labels = []
    all_pcs_2d = []
    
    for ct in common_celltypes:
        pcs = pca_projections_sample1[ct][:, :2]
        all_pcs_2d.append(pcs)
        all_labels.extend([f"{ct}_{args.sample1}"] * len(pcs))
        
        pcs = pca_projections_sample2[ct][:, :2]
        all_pcs_2d.append(pcs)
        all_labels.extend([f"{ct}_{args.sample2}"] * len(pcs))
    
    all_pcs_2d = np.vstack(all_pcs_2d)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    unique_labels = list(set(all_labels))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))
    
    for i, label in enumerate(unique_labels):
        mask = np.array(all_labels) == label
        ax.scatter(all_pcs_2d[mask, 0], all_pcs_2d[mask, 1], 
                  c=[colors[i]], label=label, alpha=0.6, s=10)
    
    ax.set_xlabel('PC1', fontsize=10)
    ax.set_ylabel('PC2', fontsize=10)
    ax.set_title(f'PCA Projection - All Cell Types\n{args.prefix}', fontsize=12, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    pca_file = os.path.join("figures", f"{args.prefix}_pca_all_celltypes.png")
    plt.savefig(pca_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {pca_file}")
    plt.close()
    
    # 3. Diagonal similarity bar plot
    fig, ax = plt.subplots(figsize=(max(8, n_ct * 0.5), 6))
    diagonal_sims = [similarity_matrix[i, i] for i in range(n_ct)]
    bars = ax.bar(common_celltypes, diagonal_sims, color=plt.cm.viridis(np.linspace(0.2, 0.8, n_ct)))
    
    ax.set_xlabel('Cell Type', fontsize=10)
    ax.set_ylabel('Cross-Sample Similarity', fontsize=10)
    ax.set_title(f'Same Cell Type Across Samples\n{args.sample1} vs {args.sample2}', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1.05])
    ax.tick_params(axis='x', rotation=45, labelsize=9)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
    ax.grid(alpha=0.3, axis='y')
    
    for bar, val in zip(bars, diagonal_sims):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
               f'{val:.3f}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    bar_file = os.path.join("figures", f"{args.prefix}_diagonal_similarities.png")
    plt.savefig(bar_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {bar_file}")
    plt.close()
    
    print("\n" + "="*50)
    print(f"Same Cell Type Similarity: {args.sample1} vs {args.sample2}")
    print("="*50)
    for i, ct in enumerate(common_celltypes):
        print(f"{ct:20s}: {similarity_matrix[i, i]:.4f}")
    print("="*50)
    print(f"\nAverage: {np.mean(diagonal_sims):.4f} ± {np.std(diagonal_sims):.4f}")

if __name__ == "__main__":
    main()
