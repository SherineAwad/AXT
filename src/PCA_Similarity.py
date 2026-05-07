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
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCA components")
    parser.add_argument("--max_cells", type=int, default=10000, help="Max cells per cell type")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if "celltype" not in adata.obs:
        raise ValueError("Missing celltype column")

    celltypes = adata.obs["celltype"].unique()
    print(f"Found cell types: {celltypes}")

    # Store PCA projections per cell type
    pca_projections = {}
    
    # Also store all cells combined for global PCA (optional)
    all_cells_by_ct = {}

    for ct in celltypes:
        print(f"Processing: {ct}")
        
        # Subset to this cell type
        idx = adata.obs["celltype"] == ct
        sub = adata[idx].copy()
        
        # Subsample if too many cells
        if sub.n_obs > args.max_cells:
            rng = np.random.default_rng(42)
            sample_idx = rng.choice(sub.n_obs, args.max_cells, replace=False)
            sub = sub[sample_idx].copy()
        
        # Get expression matrix
        X = get_X(sub)
        
        # Remove zero-variance genes
        var = np.var(X, axis=0)
        keep = var > 1e-8
        X = X[:, keep]
        
        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # PCA
        pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
        X_pca = pca.fit_transform(X_scaled)
        
        # Store
        pca_projections[ct] = X_pca
        all_cells_by_ct[ct] = X_pca
        
        print(f"  {sub.n_obs} cells, {X_pca.shape[1]} PCs kept")

    # Calculate pairwise Wasserstein distances between all cell types
    celltype_list = list(pca_projections.keys())
    n_ct = len(celltype_list)
    distance_matrix = np.zeros((n_ct, n_ct))
    
    print("\nCalculating pairwise distances...")
    
    for i, ct1 in enumerate(celltype_list):
        for j, ct2 in enumerate(celltype_list):
            if i <= j:
                # If same cell type, distance is 0
                if i == j:
                    distance_matrix[i, j] = 0
                else:
                    # Calculate Wasserstein distance for each PC dimension, then average
                    proj1 = pca_projections[ct1]
                    proj2 = pca_projections[ct2]
                    
                    n_pcs = min(proj1.shape[1], proj2.shape[1])
                    
                    # Average Wasserstein distance across PC dimensions
                    distances = []
                    for pc in range(n_pcs):
                        d = wasserstein_distance(proj1[:, pc], proj2[:, pc])
                        distances.append(d)
                    
                    distance_matrix[i, j] = np.mean(distances)
                    distance_matrix[j, i] = distance_matrix[i, j]
                    
        print(f"  Done: {ct1}")

    # Normalize distances to 0-1 range for interpretability
    max_dist = distance_matrix.max()
    if max_dist > 0:
        similarity_matrix = 1 - (distance_matrix / max_dist)
    else:
        similarity_matrix = np.ones_like(distance_matrix)

    # -------------------------
    # PLOTS
    # -------------------------
    os.makedirs("figures", exist_ok=True)
    
    # 1. Heatmap of similarity (1 - normalized distance) - CHANGED TO RED-WHITE-BLUE
    fig, ax = plt.subplots(figsize=(max(8, n_ct * 0.6), max(6, n_ct * 0.5)))
    
    im = ax.imshow(similarity_matrix, cmap='RdBu_r', aspect='auto', vmin=0, vmax=1)
    
    ax.set_xticks(range(n_ct))
    ax.set_yticks(range(n_ct))
    ax.set_xticklabels(celltype_list, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(celltype_list, fontsize=9)
    
    for i in range(n_ct):
        for j in range(n_ct):
            text = ax.text(j, i, f'{similarity_matrix[i, j]:.2f}',
                          ha="center", va="center", color="white" if similarity_matrix[i, j] < 0.5 else "black",
                          fontsize=8)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Similarity Score (1 = identical, 0 = completely different)', fontsize=10)
    plt.title(f'Cell Type Transcriptional Similarity\n{args.prefix}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    sim_heatmap = os.path.join("figures", f"{args.prefix}_celltype_similarity_heatmap.png")
    plt.savefig(sim_heatmap, dpi=300, bbox_inches='tight')
    print(f"Saved: {sim_heatmap}")
    plt.close()
    
    # 2. PCA visualization of all cell types
    all_labels = []
    all_pcs_2d = []
    
    for ct in celltype_list:
        # Take first 2 PCs
        pcs = pca_projections[ct][:, :2]
        all_pcs_2d.append(pcs)
        all_labels.extend([ct] * len(pcs))
    
    all_pcs_2d = np.vstack(all_pcs_2d)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    unique_ct = list(celltype_list)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_ct)))
    
    for i, ct in enumerate(unique_ct):
        mask = np.array(all_labels) == ct
        ax.scatter(all_pcs_2d[mask, 0], all_pcs_2d[mask, 1], 
                  c=[colors[i]], label=ct, alpha=0.6, s=10)
    
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

    print(f"\nDone: {args.output}")

if __name__ == "__main__":
    main()
