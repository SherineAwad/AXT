import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import palantir
from palantir.core import run_palantir
from palantir.presults import compute_gene_trends, cluster_gene_trends

if __name__ == '__main__':
    # -------------------------
    # args
    # -------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster_key", default="celltype", help="Column name in adata.obs for clustering")
    parser.add_argument("--root", required=True, help="Cell type to use as root (must match a value in cluster_key column)")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    # -------------------------
    # load data
    # -------------------------
    adata = sc.read_h5ad(args.input)
    print(f"[DEBUG] Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # -------------------------
    # use existing log1p layer for PCA
    # -------------------------
    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"]
        print("[DEBUG] Using existing log1p layer")
    else:
        print("[DEBUG] Warning: log1p layer not found, using adata.X as is")

    # -------------------------
    # ensure PCA exists
    # -------------------------
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')
        print("[DEBUG] PCA computed")
    else:
        print("[DEBUG] PCA already exists")

    # -------------------------
    # neighbors for visualization
    # -------------------------
    sc.pp.neighbors(adata)
    
    # Diffusion map with enough components for Palantir
    sc.tl.diffmap(adata, n_comps=min(50, adata.n_obs - 1))
    print("[DEBUG] Neighbors and diffusion map computed")
    
    # Palantir expects this specific key name
    adata.obsm['DM_EigenVectors_multiscaled'] = adata.obsm['X_diffmap'].copy()

    # -------------------------
    # ensure UMAP exists for visualization
    # -------------------------
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)
        print("[DEBUG] UMAP computed")

    # -------------------------
    # Palantir trajectory inference
    # -------------------------
    print(f"[DEBUG] Running Palantir with root: {args.root}...")
    
    # Use the provided root cell type
    root_cells = adata.obs[adata.obs[args.cluster_key] == args.root].index
    if len(root_cells) == 0:
        raise ValueError(f"Root cell type '{args.root}' not found in {args.cluster_key} column. Available: {adata.obs[args.cluster_key].unique()}")
    
    early_cell = root_cells[0]
    print(f"[DEBUG] Using early cell: {early_cell} (cell type: {args.root})")
    
    # Run Palantir
    pr_res = run_palantir(
        adata,
        early_cell=early_cell,
        num_waypoints=min(1200, adata.n_obs // 10),
        max_iterations=25,
        seed=42
    )
    
    # Add results to adata
    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obs['palantir_entropy'] = pr_res.entropy
    
    # Add branch probabilities
    for branch in pr_res.branch_probs.columns:
        adata.obs[f'palantir_branch_prob_{branch}'] = pr_res.branch_probs[branch].values
    
    # Add terminal annotation (pseudotime > 0.9 = terminal)
    adata.obs['terminal_annotation'] = 'non-terminal'
    terminal_cells = adata.obs['palantir_pseudotime'] > 0.9
    adata.obs.loc[terminal_cells, 'terminal_annotation'] = 'terminal'
    
    print("[DEBUG] Palantir completed")

    # -------------------------
    # generate Palantir plots with terminal annotation
    # -------------------------
    print("[DEBUG] Saving Palantir pseudotime UMAP...")
    sc.pl.umap(adata, color="palantir_pseudotime", 
               save=f"_{args.prefix}_palantir_pseudotime.png", show=False)
    
    print("[DEBUG] Saving terminal annotation on UMAP...")
    sc.pl.umap(adata, color="terminal_annotation", 
               save=f"_{args.prefix}_terminal_on_umap.png", show=False)
    
    print("[DEBUG] Saving Palantir entropy UMAP...")
    sc.pl.umap(adata, color="palantir_entropy",
               save=f"_{args.prefix}_palantir_entropy.png", show=False)
    
    # -------------------------
    # Plot branch probabilities for each terminal branch
    # -------------------------
    branch_cols = [c for c in adata.obs.columns if c.startswith("palantir_branch_prob_")]
    if len(branch_cols) > 0:
        print(f"[DEBUG] Found {len(branch_cols)} terminal branches")
        for branch_col in branch_cols:
            branch_name = branch_col.replace("palantir_branch_prob_", "")
            sc.pl.umap(adata, color=branch_col, cmap="viridis",
                       title=f"Probability: {branch_name}",
                       save=f"_{args.prefix}_branch_{branch_name}.png", show=False)
    
    # -------------------------
    # PAGA plots with pseudotime and terminal annotation
    # -------------------------
    if args.cluster_key in adata.obs.columns:
        print("[DEBUG] Computing PAGA...")
        sc.tl.paga(adata, groups=args.cluster_key)
        
        print("[DEBUG] Saving PAGA with pseudotime...")
        sc.pl.paga(adata, color="palantir_pseudotime", 
                   save=f"_{args.prefix}_paga_pseudotime.png", show=False)
        
        print("[DEBUG] Saving PAGA with terminal annotation...")
        sc.pl.paga(adata, color="terminal_annotation", 
                   save=f"_{args.prefix}_paga_terminal.png", show=False)

    # -------------------------
    # Print terminal branch information
    # -------------------------
    print("\n" + "="*50)
    print("PALANTIR RESULTS")
    print("="*50)
    
    # Identify terminal branches
    branch_cols = [c for c in adata.obs.columns if c.startswith("palantir_branch_prob_")]
    if len(branch_cols) > 0:
        print(f"Number of terminal branches detected: {len(branch_cols)}")
        print("\nTerminal branches and their dominant cell types:")
        
        for branch_col in branch_cols:
            branch_name = branch_col.replace("palantir_branch_prob_", "")
            # Find cells with high probability (>0.8) for this branch
            high_prob_cells = adata.obs[adata.obs[branch_col] > 0.8]
            if len(high_prob_cells) > 0:
                # Get most common cell type in this branch
                celltype_counts = high_prob_cells[args.cluster_key].value_counts()
                dominant_celltype = celltype_counts.index[0]
                print(f"  Branch '{branch_name}': {dominant_celltype} ({len(high_prob_cells)} cells with >0.8 probability)")
            else:
                print(f"  Branch '{branch_name}': No cells >0.8 probability")
    else:
        print("No branch probabilities found - linear trajectory detected")
    
    # Print pseudotime range
    if "palantir_pseudotime" in adata.obs.columns:
        pt_min = adata.obs["palantir_pseudotime"].min()
        pt_max = adata.obs["palantir_pseudotime"].max()
        print(f"\nPseudotime range: {pt_min:.2f} to {pt_max:.2f}")
    
    # Count terminal vs non-terminal cells
    n_terminal = terminal_cells.sum()
    n_non_terminal = (~terminal_cells).sum()
    print(f"\nTerminal cells (pseudotime > 0.9): {n_terminal}")
    print(f"Non-terminal cells: {n_non_terminal}")
    
    print("="*50 + "\n")

    # -------------------------
    # summary
    # -------------------------
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"Total cells analyzed: {adata.n_obs}")
    print(f"Root cell type: {args.root}")
    print(f"Number of terminal branches: {len(branch_cols)}")
    print(f"\nOutputs saved to figures/ directory:")
    print(f"  - Palantir pseudotime: figures/umap_{args.prefix}_palantir_pseudotime.png")
    print(f"  - Terminal annotation on UMAP: figures/umap_{args.prefix}_terminal_on_umap.png")
    print(f"  - Palantir entropy: figures/umap_{args.prefix}_palantir_entropy.png")
    print(f"  - PAGA with pseudotime: figures/paga_{args.prefix}_paga_pseudotime.png")
    print(f"  - PAGA with terminal annotation: figures/paga_{args.prefix}_paga_terminal.png")
    print(f"  - Branch probability UMAPs per terminal branch")
    print(f"\nOutput file: {args.output}")
    print("="*50)
