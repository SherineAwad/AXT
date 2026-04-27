#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--remove", required=False)
parser.add_argument("--min_exp", type=float, default=0, required=False)
parser.add_argument("--max_egfp", type=float, default=2.0, required=False)
args = parser.parse_args()

# save folder
sc.settings.figdir = "figures"

# load
adata = sc.read(args.input)

# -------------------------
# REMOVE MASK
# -------------------------
remove_mask = np.zeros(adata.n_obs, dtype=bool)

if args.remove:
    # Cdkn2a expression
    gene = args.remove
    gene_idx = adata.var_names.get_loc(gene)
    
    if hasattr(adata.X, "toarray"):
        expr_cdkn2a = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr_cdkn2a = adata.X[:, gene_idx]
    
    # EGFP expression
    egfp_idx = adata.var_names.get_loc("EGFP")
    
    if hasattr(adata.X, "toarray"):
        expr_egfp = adata.X[:, egfp_idx].toarray().flatten()
    else:
        expr_egfp = adata.X[:, egfp_idx]
    
    cell_mask = adata.obs["celltype"] == args.celltype
    
    # Condition: Cdkn2a > min_exp AND EGFP < max_egfp
    remove_mask = cell_mask & (expr_cdkn2a > args.min_exp) & (expr_egfp < args.max_egfp)

# Save original counts
original_count = adata.n_obs
original_celltype_count = (adata.obs["celltype"] == args.celltype).sum()

# -------------------------
# BEFORE UMAP - show what will be removed
# -------------------------
adata.obs["temp_remove"] = "keep"
adata.obs.loc[remove_mask, "temp_remove"] = "will remove"

# Create filename with parameters
fig_before = f"_{args.prefix}_before_minExp{args.min_exp}_maxEGFP{args.max_egfp}.png"

sc.pl.umap(
    adata,
    color="temp_remove",
    show=False,
    save=fig_before
)

# -------------------------
# APPLY FILTER
# -------------------------
adata = adata[~remove_mask, :].copy()

# After counts
after_count = adata.n_obs
after_celltype_count = (adata.obs["celltype"] == args.celltype).sum()

# -------------------------
# VERIFICATION
# -------------------------
print(f"Total cells before: {original_count}")
print(f"Total cells after: {after_count}")
print(f"Total cells removed: {original_count - after_count}")
print(f"")
print(f"{args.celltype} cells before: {original_celltype_count}")
print(f"{args.celltype} cells after: {after_celltype_count}")
print(f"{args.celltype} cells removed: {original_celltype_count - after_celltype_count}")

# -------------------------
# AFTER UMAP
# -------------------------
fig_after = f"_{args.prefix}_after_minExp{args.min_exp}_maxEGFP{args.max_egfp}.png"

sc.pl.umap(
    adata,
    color="celltype",
    show=False,
    save=fig_after
)

# -------------------------
# SAVE OBJECT (ALWAYS LAST)
# -------------------------
adata.write(args.output)

