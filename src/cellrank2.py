#!/usr/bin/env python3

import argparse
import scanpy as sc
import cellrank as cr
import matplotlib.pyplot as plt
import scvelo as scv
import palantir
import numpy as np
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
    
    adata.obs['macrostates'] = pt_est.macrostates
    
    palantir.plot.plot_trajectories(adata, cell_color='macrostates')
    plt.title(f"Macrostates - Root: {root_name}")
    plt.savefig(f"figures/{args.prefix}_{root_name}_macrostates.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    pt_est.predict_terminal_states()
    adata.obs['term_states_fwd'] = pt_est.terminal_states
    
    terminal_names = adata.obs['term_states_fwd'].cat.categories.tolist()
    
    palantir.plot.plot_trajectories(adata, cell_color='term_states_fwd')
    plt.title(f"Terminal States - Root: {root_name}")
    plt.figtext(0.5, 0.01, f"Terminal states: {', '.join(terminal_names)}", ha="center", fontsize=8, bbox=dict(facecolor='white', alpha=0.8))
    plt.savefig(f"figures/{args.prefix}_{root_name}_terminal_states.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    pt_est.compute_fate_probabilities()
    
    for lineage in pt_est.fate_probabilities.names:
        adata.obs[f'fate_prob_{root_name}_{lineage}'] = np.array(pt_est.fate_probabilities[lineage])
    
    # ONE FIGURE with all fate probabilities as panels
    lineages = pt_est.fate_probabilities.names
    n_lineages = len(lineages)
    ncols = 3
    nrows = (n_lineages + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    if n_lineages > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    
    for i, lineage in enumerate(lineages):
        scv.pl.scatter(adata, color=f'fate_prob_{root_name}_{lineage}', 
                      basis='X_umap', ax=axes[i], show=False, title=lineage)
    
    for i in range(n_lineages, len(axes)):
        axes[i].axis('off')
    
    plt.suptitle(f"Fate Probabilities - Root: {root_name}", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"figures/{args.prefix}_{root_name}_fate_probabilities.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    # Driver genes (one figure per lineage, multiple gene panels)
    for lineage in pt_est.fate_probabilities.names:
        drivers = pt_est.compute_lineage_drivers(lineages=[lineage], cluster_key='celltype')
        drivers_sorted = drivers.sort_values(by=f'{lineage}_corr', ascending=False)
        top_genes = drivers_sorted.index.tolist()[1:10]
        sc.pl.umap(adata, color=top_genes, legend_loc='on data', legend_fontsize='x-small', use_raw=False, ncols=3)
        plt.savefig(f"figures/{args.prefix}_{root_name}_drivers_{lineage}.png", dpi=150, bbox_inches="tight")
        plt.close()

adata.write_h5ad(args.output)
print(f"Done. Output saved to {args.output}")
