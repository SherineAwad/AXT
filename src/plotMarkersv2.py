#!/usr/bin/env python3

import argparse
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--markers", required=True)
args = parser.parse_args()

# Read marker file
marker_genes = {}

with open(args.markers) as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        celltype, genes = line.split(":")
        marker_genes[celltype] = [g.strip() for g in genes.split(",")]

adata = sc.read_h5ad(args.input)

# Collect all unique genes
all_genes = []
for genes in marker_genes.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))

# Filter and warn about missing genes
valid_genes = []
for gene in all_genes:
    if gene in adata.var_names:
        valid_genes.append(gene)
    else:
        print(f"WARNING: {gene} not found in dataset - skipping")

# Feature plot for each gene (gene name only in filename)
for gene in valid_genes:
    sc.pl.umap(adata, color=gene, save=f"_{args.prefix}_{gene}.png")

# Dotplot (uses full marker_genes dict with celltypes)
filtered_marker_genes = {}
for celltype, genes in marker_genes.items():
    filtered_genes = [g for g in genes if g in adata.var_names]
    if filtered_genes:
        filtered_marker_genes[celltype] = filtered_genes

sc.pl.dotplot(
    adata,
    filtered_marker_genes,
    groupby="leiden",
    save=f"_{args.prefix}_dotplot.png"
)


