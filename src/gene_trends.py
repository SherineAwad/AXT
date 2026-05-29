import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import gc
import palantir
from palantir.presults import compute_gene_trends

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Saved Palantir h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file (with gene trends added)")
    parser.add_argument("--prefix", required=True, help="Prefix for output plot files")
    parser.add_argument("--cluster_key", default="leiden", help="Column name in adata.obs for clustering")
    parser.add_argument("--root", required=False, default=None, help="Specific root value")
    parser.add_argument("--n_top_genes", type=int, default=20, help="Number of top varying genes to plot")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    adata = sc.read_h5ad(args.input)
    print(f"[DEBUG] Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Determine which roots to process
    if args.root is not None and args.root != "":
        roots_to_use = [args.root]
        print(f"[DEBUG] Single root mode: {args.root}")
    else:
        pseudotime_cols = [c for c in adata.obs.columns if c.startswith("palantir_pseudotime_")]
        if len(pseudotime_cols) > 0:
            roots_to_use = [c.replace("palantir_pseudotime_", "") for c in pseudotime_cols]
        else:
            roots_to_use = adata.obs[args.cluster_key].unique()
            roots_to_use = [r for r in roots_to_use if not pd.isna(r)]
        print(f"[DEBUG] Processing {len(roots_to_use)} roots: {roots_to_use}")

    for current_root in roots_to_use:
        print("\n" + "="*50)
        print(f"[DEBUG] Processing root: {current_root}")

        pseudotime_col = f"palantir_pseudotime_{current_root}"
        if pseudotime_col not in adata.obs.columns:
            if "palantir_pseudotime" in adata.obs.columns and len(roots_to_use) == 1:
                pseudotime_col = "palantir_pseudotime"
            else:
                print(f"[WARNING] Pseudotime column '{pseudotime_col}' not found. Skipping.")
                continue

        print(f"[DEBUG] Using pseudotime column: {pseudotime_col}")

        # Store original pseudotime, then rename for compute_gene_trends
        original_pseudotime = adata.obs[pseudotime_col].copy()
        adata.obs['palantir_pseudotime'] = original_pseudotime

        valid_cells = ~pd.isna(adata.obs['palantir_pseudotime'])
        adata_subset = adata[valid_cells].copy()
        print(f"[DEBUG] Using {adata_subset.n_obs} cells with valid pseudotime")

        try:
            print("[DEBUG] Computing gene expression trends...")
            gene_trends_all = compute_gene_trends(adata_subset, expression_key=None)
            
            # Convert to DataFrame if dict
            if isinstance(gene_trends_all, dict):
                print("[DEBUG] Converting dictionary to DataFrame...")
                if 'trends' in gene_trends_all:
                    gene_trends_df = gene_trends_all['trends']
                else:
                    gene_trends_df = pd.DataFrame.from_dict(gene_trends_all, orient='index')
            else:
                gene_trends_df = gene_trends_all
            
            print(f"[DEBUG] Gene trends computed for {len(gene_trends_df)} genes")
            
            # Store full trends in adata.uns
            adata.uns[f'gene_trends_root_{current_root}'] = gene_trends_df

            # Get top varying genes by variance
            trend_vars = gene_trends_df.var(axis=1).sort_values(ascending=False)
            top_genes = trend_vars.head(args.n_top_genes).index.tolist()
            print(f"[DEBUG] Top {len(top_genes)} genes for root {current_root}: {top_genes[:5]}...")

            if len(top_genes) == 0:
                print("[WARNING] No genes found for plotting. Skipping plot.")
                continue

            # Plot gene trends (without 'ax' parameter for older Palantir version)
            print("[DEBUG] Plotting gene trends...")
            plt.figure(figsize=(14, 10))
            palantir.plot.plot_gene_trends(adata_subset, genes=top_genes)
            plt.title(f"Gene expression trends (root: {current_root})")
            plt.tight_layout()
            plt.savefig(f"figures/{args.prefix}_root_{current_root}_gene_trends.png", dpi=150, bbox_inches="tight")
            plt.close()
            print(f"[DEBUG] Saved: figures/{args.prefix}_root_{current_root}_gene_trends.png")

            # Free memory
            del gene_trends_all, gene_trends_df
            gc.collect()

        except Exception as e:
            print(f"[WARNING] Could not compute or plot gene trends for root {current_root}: {e}")

    # Save updated h5ad
    adata.write_h5ad(args.output)
    print(f"\n[DEBUG] Saved annotated data to {args.output}")

    print("\n" + "="*50)
    print("GENE TREND ANALYSIS COMPLETED")
    print("="*50)
    print(f"Processed {len(roots_to_use)} roots")
    print(f"Output h5ad: {args.output}")
    print("="*50)
