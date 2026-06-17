#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action="ignore", category=Warning)

sc.settings.verbosity = 3
sc.settings.set_figure_params(
    dpi=100,
    frameon=False,
    figsize=(5, 5),
    facecolor='white',
    color_map='viridis_r'
)

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--sample", required=False, default=None)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# OPTIONAL SAMPLE FILTER
suffix = ""
if args.sample is not None:
    adata = adata[adata.obs["sample"] == args.sample].copy()
    suffix = f"_sample-{args.sample}"


# Plot : PAGA abstract cluster graph
sc.tl.paga(adata, groups='celltype')
sc.pl.paga(
    adata,
    color='celltype',
    edge_width_scale=0.3,
    show=False
)
plt.savefig(
    f"figures/{args.prefix}{suffix}_paga_clusters.png",
    bbox_inches='tight',
    dpi=150
)
plt.close()

print(f"All plots saved to figures/{args.prefix}{suffix}_*.png")
