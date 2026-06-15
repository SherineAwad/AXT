import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import palantir
from palantir.core import run_palantir

if __name__ == '__main__':
    # -------------------------
    # args
    # -------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster_key", default="celltype", help="Column name in adata.obs for clustering")
    parser.add_argument("--root", required=False, default=None, help="Specific root value from cluster_key. If not provided, loops through all unique values.")
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
    # Determine which roots to use
    # -------------------------
    if args.root is not None and args.root != "":
        roots_to_use = [args.root]
        print(f"[DEBUG] Single root mode: {args.root}")
    else:
        roots_to_use = adata.obs[args.cluster_key].unique()
        roots_to_use = [r for r in roots_to_use if not pd.isna(r)]
        print(f"[DEBUG] Loop mode: Running Palantir for each of {len(roots_to_use)} roots: {roots_to_use}")

    # -------------------------
    # Store results for all roots to create summary table
    # -------------------------
    all_results_summary = []

    # -------------------------
    # Run Palantir for each root
    # -------------------------
    for current_root in roots_to_use:
        print("\n" + "="*50)
        print(f"[DEBUG] Running Palantir with root: {current_root}")
        print("="*50)
        
        root_cells = adata.obs[adata.obs[args.cluster_key] == current_root].index
        if len(root_cells) == 0:
            print(f"[WARNING] Root '{current_root}' not found in {args.cluster_key} column. Skipping.")
            continue
        
        early_cell = root_cells[0]
        print(f"[DEBUG] Using early cell: {early_cell} (root: {current_root})")
        
        # Run Palantir
        pr_res = run_palantir(
            adata,
            early_cell=early_cell,
            num_waypoints=min(1200, adata.n_obs // 10),
            max_iterations=25,
            seed=42
        )
        
        # Add results to adata with root-specific column names
        adata.obs[f'palantir_pseudotime_{current_root}'] = pr_res.pseudotime
        adata.obs[f'palantir_entropy_{current_root}'] = pr_res.entropy
        
        # Add branch probabilities with root-specific names
        for branch in pr_res.branch_probs.columns:
            adata.obs[f'palantir_branch_prob_{current_root}_{branch}'] = pr_res.branch_probs[branch].values
        
        # Normalize branch probabilities to sum to 1
        branch_cols = [c for c in adata.obs.columns if c.startswith(f"palantir_branch_prob_{current_root}_")]
        if len(branch_cols) > 0:
            prob_sum = adata.obs[branch_cols].sum(axis=1)
            prob_sum = prob_sum.replace(0, 1)
            for col in branch_cols:
                adata.obs[col] = adata.obs[col] / prob_sum
            print("[DEBUG] Branch probabilities normalized to sum to 1")
        
        # Select branch cells for trajectory plotting (from tutorial)
        try:
            # FIX: pass pr_res as first argument, and use the returned masks
            masks = palantir.presults.select_branch_cells(pr_res, adata, q=0.01, eps=0.01)
            print("[DEBUG] Branch masks computed using select_branch_cells")
        except Exception as e:
            print(f"[WARNING] Could not compute branch masks: {e}")
            masks = None
        
        print("[DEBUG] Palantir completed")
        
        root_prefix = f"{args.prefix}_{current_root}"
        
        # -------------------------
        # Plot: Trajectories using branch masks (from tutorial)
        # -------------------------
        try:
            print("[DEBUG] Plotting trajectories with branch masks...")
            fig, ax = plt.subplots(figsize=(10, 8))
            # Use the root-specific pseudotime for coloring
            adata.obs['temp_pseudotime'] = adata.obs[f'palantir_pseudotime_{current_root}']
            # FIX: pass branch_masks=masks to actually draw trajectories
            palantir.plot.plot_trajectories(adata, branch_masks=masks, cell_color="temp_pseudotime")
            plt.title(f"Trajectories (root: {current_root})")
            plt.savefig(f"figures/{root_prefix}_trajectories.png", dpi=150, bbox_inches="tight")
            plt.close()
            # Clean up temporary column
            del adata.obs['temp_pseudotime']
            print(f"[DEBUG] Trajectories plot saved: figures/{root_prefix}_trajectories.png")
        except Exception as e:
            print(f"[WARNING] Could not plot trajectories: {e}")
        
        # -------------------------
        # Calculate entropy based on cell assignments
        # -------------------------
        branch_cols = [c for c in adata.obs.columns if c.startswith(f"palantir_branch_prob_{current_root}_")]
        branch_assignment_entropy = 0.0
        
        if len(branch_cols) > 0:
            # Get cell assignments (which branch each cell belongs to)
            branch_prob_matrix = adata.obs[branch_cols].values
            branch_assignments = np.argmax(branch_prob_matrix, axis=1)
            
            # Count cells per branch
            cells_per_branch = []
            for i in range(len(branch_cols)):
                count = (branch_assignments == i).sum()
                cells_per_branch.append(count)
            
            # Calculate entropy from cell counts
            total_cells = sum(cells_per_branch)
            branch_probs_from_counts = [c/total_cells for c in cells_per_branch if c > 0]
            
            if len(branch_probs_from_counts) > 1:
                branch_assignment_entropy = -sum(p * np.log(p) for p in branch_probs_from_counts)
            else:
                branch_assignment_entropy = 0.0
            
            print(f"[DEBUG] Cells per branch: {cells_per_branch}")
            print(f"[DEBUG] Entropy from cell assignments: {branch_assignment_entropy:.6f}")
        
        # -------------------------
        # Print and save terminal branch information
        # -------------------------
        print("\n" + "="*50)
        print("PALANTIR RESULTS")
        print("="*50)
        
        results_list = []
        branch_cols = [c for c in adata.obs.columns if c.startswith(f"palantir_branch_prob_{current_root}_")]
        
        if len(branch_cols) > 0:
            print(f"Number of terminal branches detected: {len(branch_cols)}")
            results_list.append({"Metric": "Number of terminal branches", "Value": len(branch_cols)})
            
            branch_prob_matrix = adata.obs[branch_cols].values
            row_sums = branch_prob_matrix.sum(axis=1)
            print(f"Branch probability sum check - min: {row_sums.min():.3f}, max: {row_sums.max():.3f}")
            results_list.append({"Metric": "Branch probability sum min", "Value": f"{row_sums.min():.3f}"})
            results_list.append({"Metric": "Branch probability sum max", "Value": f"{row_sums.max():.3f}"})
            
            print("\nTerminal branches:")
            for i, branch_col in enumerate(branch_cols):
                branch_name = branch_col.replace(f"palantir_branch_prob_{current_root}_", "")
                branch_assignments = np.argmax(branch_prob_matrix, axis=1)
                cells_in_branch = (branch_assignments == i).sum()
                
                branch_cells = adata.obs[branch_assignments == i]
                if len(branch_cells) > 0:
                    celltype_counts = branch_cells[args.cluster_key].value_counts()
                    dominant_celltype = celltype_counts.index[0] if len(celltype_counts) > 0 else "none"
                    print(f"  Branch '{branch_name}': {dominant_celltype} ({cells_in_branch} cells assigned)")
                    results_list.append({"Metric": f"Branch '{branch_name}' dominant cell type", "Value": dominant_celltype})
                    results_list.append({"Metric": f"Branch '{branch_name}' cells assigned", "Value": cells_in_branch})
                else:
                    print(f"  Branch '{branch_name}': no cells assigned")
                    results_list.append({"Metric": f"Branch '{branch_name}' dominant cell type", "Value": "none"})
                    results_list.append({"Metric": f"Branch '{branch_name}' cells assigned", "Value": 0})
        else:
            print("No branch probabilities found - linear trajectory detected")
            results_list.append({"Metric": "Trajectory type", "Value": "linear"})
        
        pt_col = f'palantir_pseudotime_{current_root}'
        if pt_col in adata.obs.columns:
            pt_min = adata.obs[pt_col].min()
            pt_max = adata.obs[pt_col].max()
            print(f"\nPseudotime range: {pt_min:.2f} to {pt_max:.2f}")
            results_list.append({"Metric": "Pseudotime min", "Value": pt_min})
            results_list.append({"Metric": "Pseudotime max", "Value": pt_max})
        
        # Add entropy from cell assignments
        results_list.append({"Metric": "Branch entropy (from cell counts)", "Value": branch_assignment_entropy})
        results_list.append({"Metric": "Branch distribution (cells per branch)", "Value": str([int(x) for x in cells_per_branch]) if len(branch_cols) > 0 else "linear"})
        
        results_list.append({"Metric": "Root", "Value": str(current_root)})
        results_list.append({"Metric": "Total cells", "Value": adata.n_obs})
        
        print("="*50 + "\n")
        
        # Save results to CSV for this root
        results_df = pd.DataFrame(results_list)
        csv_filename = f"{root_prefix}_palantir_summary.csv"
        results_df.to_csv(csv_filename, index=False)
        print(f"[DEBUG] Summary saved to: {csv_filename}")
        
        # Store summary for final comparison
        all_results_summary.append({
            "Root": str(current_root),
            "Num_branches": len(branch_cols),
            "Branch_entropy_from_counts": branch_assignment_entropy,
            "Cells_per_branch": str([int(x) for x in cells_per_branch]) if len(branch_cols) > 0 else "linear",
            "Pseudotime_min": pt_min if pt_col in adata.obs.columns else np.nan,
            "Pseudotime_max": pt_max if pt_col in adata.obs.columns else np.nan
        })
        
        # Summary for this root
        print("\n" + "="*50)
        print(f"SUMMARY FOR ROOT: {current_root}")
        print("="*50)
        print(f"Total cells analyzed: {adata.n_obs}")
        print(f"Root: {current_root}")
        print(f"Number of terminal branches: {len(branch_cols)}")
        print(f"Branch entropy (from cell counts): {branch_assignment_entropy:.6f}")
        if len(branch_cols) > 0:
            print(f"Cells per branch: {[int(x) for x in cells_per_branch]}")
        print(f"\nOutputs saved:")
        print(f"  - Cluster UMAP: figures/umap_{args.prefix}_{args.cluster_key}.png")
        print(f"  - Trajectories plot: figures/{root_prefix}_trajectories.png")
        print(f"  - CSV summary: {csv_filename}")
        print("="*50)
    
    # Save annotated data with all root-specific columns
    adata.write_h5ad(args.output)
    print(f"[DEBUG] Saved annotated data to {args.output}")
    
    # -------------------------
    # Create summary table comparing all roots
    # -------------------------
    if len(all_results_summary) > 1:
        print("\n" + "="*70)
        print("SUMMARY TABLE: COMPARING ALL ROOTS")
        print("="*70)
        summary_df = pd.DataFrame(all_results_summary)
        print(summary_df.to_string(index=False))
        
        # Save summary table to CSV
        summary_csv = f"{args.prefix}_all_roots_comparison.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"\n[DEBUG] Summary table saved to: {summary_csv}")
        
        # Determine most probable trajectory based on entropy from cell counts
        if "Branch_entropy_from_counts" in summary_df.columns:
            branching_roots = summary_df[summary_df["Branch_entropy_from_counts"] > 0.1]
            if len(branching_roots) > 0:
                best_branching_root = branching_roots.loc[branching_roots["Branch_entropy_from_counts"].idxmax(), "Root"]
                best_entropy = branching_roots["Branch_entropy_from_counts"].max()
                print(f"\n[INFO] Most probable branching trajectory: Root {best_branching_root} (entropy: {best_entropy:.4f})")
            else:
                print("\n[INFO] No roots show evidence of branching (all entropy near 0).")
                print("       The data most likely follows a linear trajectory.")
            
            linear_roots = summary_df[summary_df["Branch_entropy_from_counts"] == 0]["Root"].tolist()
            if len(linear_roots) > 0:
                print(f"\n[INFO] Roots producing linear trajectory (1 effective branch): {linear_roots}")
    else:
        print("\n[INFO] Only one root analyzed. No comparison needed.")
    
    print("\n" + "="*50)
    print("ALL ROOTS COMPLETED")
    print("="*50)
    print(f"Processed {len(roots_to_use)} roots")
    print(f"Output file: {args.output}")
    print("="*50)
