#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action="ignore", category=Warning)

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(5, 5), facecolor='white', color_map='viridis_r')

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# Plot 1: UMAP overview
sc.pl.umap(adata, color=['celltype','leiden'], legend_loc='on data', legend_fontsize='xx-small', ncols=2)
plt.savefig(f"figures/{args.prefix}_umap_overview.png", bbox_inches='tight', dpi=150)
plt.close()

# Plot 3: PAGA abstract cluster graph
sc.tl.paga(adata, groups='celltype')
sc.pl.paga(adata, color='celltype', edge_width_scale=0.3)
plt.savefig(f"figures/{args.prefix}_paga_clusters.png", bbox_inches='tight', dpi=150)
plt.close()

print(f"All plots saved to figures/{args.prefix}_*.png")
