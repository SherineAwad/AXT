#!/usr/bin/env python3

import argparse
import scanpy as sc
import cellrank as cr
import matplotlib.pyplot as plt
import palantir
import os
import warnings
warnings.simplefilter(action="ignore", category=Warning)

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

adata = sc.read_h5ad(args.input)

pt_columns = [col for col in adata.obs.columns if col.startswith("palantir_pseudotime_")]

for pt_col in pt_columns:
    root_name = pt_col.replace("palantir_pseudotime_", "")
    
    print(f"Processing root: {root_name}")
    
    pt_kernel = cr.kernels.PseudotimeKernel(adata, time_key=pt_col)
    pt_kernel.compute_transition_matrix()
    
    pt_est = cr.estimators.GPCCA(pt_kernel)
    pt_est.fit(cluster_key='celltype')
    
    # Get macrostates
    adata.obs['macrostates'] = pt_est.macrostates
    macro_names = adata.obs['macrostates'].cat.categories.tolist()
    
    # Plot Macrostates
    palantir.plot.plot_trajectories(adata, cell_color='macrostates')
    plt.title(f"Macrostates - Root: {root_name}")
    plt.figtext(0.5, 0.01, f"Macrostates: {', '.join(macro_names)}", ha="center", fontsize=8, bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))
    plt.savefig(f"figures/{args.prefix}_{root_name}_macrostates.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # Get terminal states
    pt_est.predict_terminal_states()
    adata.obs['term_states_fwd'] = pt_est.terminal_states
    terminal_names = adata.obs['term_states_fwd'].cat.categories.tolist()
    
    # Plot Terminal States
    palantir.plot.plot_trajectories(adata, cell_color='term_states_fwd')
    plt.title(f"Terminal States - Root: {root_name}")
    plt.figtext(0.5, 0.01, f"Terminal states: {', '.join(terminal_names)}", ha="center", fontsize=8, bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))
    plt.savefig(f"figures/{args.prefix}_{root_name}_terminal_states.png", dpi=150, bbox_inches="tight")
    plt.close()

adata.write_h5ad(args.output)
print(f"Done. Output saved to {args.output}")
