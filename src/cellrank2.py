import scanpy as sc
import cellrank as cr
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # -------------------------
    # args
    # -------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster_key", default="celltype", help="Column name in adata.obs for clustering")
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
    # neighbors + diffusion
    # -------------------------
    sc.pp.neighbors(adata)
    sc.tl.diffmap(adata)
    print("[DEBUG] Neighbors and diffusion map computed")

    # -------------------------
    # CellRank 2 trajectory - auto-detect number of terminal states
    # -------------------------
    print("[DEBUG] Computing transition matrix...")
    k = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    
    print("[DEBUG] Initializing GPCCA estimator...")
    gpi = cr.estimators.GPCCA(k)
    
    print("[DEBUG] Computing Schur decomposition...")
    # Use all components to detect natural number of states
    gpi.compute_schur(n_components=min(50, adata.n_obs - 1))
    
    # Auto-detect number of macrostates using eigengap heuristic
    print("[DEBUG] Auto-detecting number of macrostates...")
    gpi.compute_macrostates(cluster_key=args.cluster_key)
    
    # Get the auto-detected number
    n_states = len(gpi.macrostates)
    print(f"[DEBUG] Auto-detected {n_states} macrostates")
    
    print("[DEBUG] Predicting terminal states...")
    gpi.predict_terminal_states()
    
    print("[DEBUG] Computing fate probabilities...")
    gpi.compute_fate_probabilities()
    
    print("[DEBUG] Computing lineage drivers...")
    gpi.compute_lineage_drivers(lineages=None, key_added="lineage_drivers")

    # -------------------------
    # ensure UMAP exists
    # -------------------------
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)
        print("[DEBUG] UMAP computed")

    # -------------------------
    # generate plots
    # -------------------------
    print("[DEBUG] Saving lineage UMAP...")
    sc.pl.umap(adata, color="macrostates_fwd", save=f"_{args.prefix}_lineages.png", show=False)

    # -------------------------
    # CellRank's own plots
    # -------------------------
    if args.cluster_key in adata.obs.columns:
        print("[DEBUG] Computing PAGA...")
        sc.tl.paga(adata, groups=args.cluster_key)
        
        print("[DEBUG] Saving aggregate fate probabilities plot...")
        cr.pl.aggregate_fate_probabilities(adata, basis="umap", cluster_key=args.cluster_key, save=f"{args.prefix}_cellrank_fates.png", show=False)
        
        print("[DEBUG] Saving PAGA graph...")
        sc.pl.paga(adata, color=args.cluster_key, save=f"_{args.prefix}_paga_graph.png", show=False)

    # -------------------------
    # Print terminal states with their cell type names
    # -------------------------
    print("\n" + "="*50)
    print("TERMINAL STATES IDENTIFIED")
    print("="*50)
    
    if hasattr(gpi, 'terminal_states') and gpi.terminal_states is not None:
        terminal_list = list(gpi.terminal_states)
        print(f"Number of terminal states: {len(terminal_list)}")
        print("\nTerminal states and their dominant cell types:")
        
        # Map each terminal macrostate to cell type
        unique_states = adata.obs['macrostates_fwd'].dropna().unique()
        for state in sorted(unique_states):
            if state in terminal_list:
                cells_in_state = adata.obs[adata.obs['macrostates_fwd'] == state]
                if len(cells_in_state) > 0:
                    # Get most common cell type in this state
                    celltype_counts = cells_in_state[args.cluster_key].value_counts()
                    dominant_celltype = celltype_counts.index[0]
                    print(f"  Terminal {state}: {dominant_celltype} ({celltype_counts[dominant_celltype]} cells)")
    
    print("="*50 + "\n")

    # -------------------------
    # Print lineage drivers information
    # -------------------------
    print("\n" + "="*50)
    print("LINEAGE DRIVERS PER TERMINAL STATE")
    print("="*50)
    
    if 'terminal_lineage_drivers' in adata.varm:
        driver_df = adata.varm['terminal_lineage_drivers']
        corr_cols = [c for c in driver_df.columns if c.endswith('_corr')]
        
        # Get terminal state names from the mapping
        for i, col in enumerate(corr_cols):
            # Try to map to cell type
            if i < len(unique_states):
                state_val = sorted([s for s in unique_states if not pd.isna(s)])[i] if i < len([s for s in unique_states if not pd.isna(s)]) else i
                cells_in_state = adata.obs[adata.obs['macrostates_fwd'] == state_val] if 'state_val' in locals() else None
                if cells_in_state is not None and len(cells_in_state) > 0:
                    celltype = cells_in_state[args.cluster_key].value_counts().index[0]
                    print(f"\n{celltype} (Terminal {state_val}):")
                else:
                    print(f"\nTerminal {i}:")
            else:
                print(f"\nTerminal {i}:")
            
            top_genes = driver_df.nlargest(10, col).index.tolist()
            print(f"  Top 5 genes: {', '.join(top_genes[:5])}")
    else:
        print("Warning: 'terminal_lineage_drivers' not found")
    
    print("="*50 + "\n")

    # -------------------------
    # Heatmaps for each terminal state
    # -------------------------
    if 'terminal_lineage_drivers' in adata.varm:
        print("[DEBUG] Generating lineage driver heatmaps...")
        adata_clean = adata[~pd.isna(adata.obs['macrostates_fwd'])].copy()
        
        if adata_clean.n_obs > 0:
            driver_df = adata_clean.varm['terminal_lineage_drivers']
            corr_cols = [c for c in driver_df.columns if c.endswith('_corr')]
            
            for i, col in enumerate(corr_cols[:10]):  # Limit to first 10 for readability
                lineage_name = col.replace('_corr', '')
                top_genes = driver_df.nlargest(20, col).index.tolist()
                sc.pl.heatmap(adata_clean, var_names=top_genes, groupby="macrostates_fwd", cmap="viridis",
                              save=f"_{args.prefix}_lineage_drivers_{lineage_name}.png", show=False)
    
    # -------------------------
    # summary
    # -------------------------
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"Total cells analyzed: {adata.n_obs}")
    print(f"Auto-detected number of macrostates: {n_states}")
    print(f"Number of terminal states: {len(terminal_list) if hasattr(gpi, 'terminal_states') and gpi.terminal_states is not None else 'Unknown'}")
    print(f"\nOutputs saved to figures/ directory")
    print(f"Output file: {args.output}")
    print("="*50)
