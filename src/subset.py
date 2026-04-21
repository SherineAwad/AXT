import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input h5ad file")
parser.add_argument("--output", required=True, help="Output h5ad file")
parser.add_argument("--celltype", required=True, help="Celltype to subset (e.g., T cells)")
parser.add_argument("--prefix", default="output", help="Prefix for saved figures")
args = parser.parse_args()

# -------------------------
# Read data
# -------------------------
adata = sc.read_h5ad(args.input)

# -------------------------
# Validation
# -------------------------
if "celltype" not in adata.obs.columns:
    raise ValueError("celltype column not found in adata.obs")

if args.celltype not in adata.obs["celltype"].unique():
    available = adata.obs["celltype"].unique()
    raise ValueError(f"Celltype '{args.celltype}' not found. Available: {available}")

if "sample" not in adata.obs.columns:
    raise ValueError("sample column not found in adata.obs")

# -------------------------
# Subset
# -------------------------
adata_subset = adata[adata.obs["celltype"] == args.celltype].copy()

print(f"Subsetted {args.celltype}: {adata_subset.n_obs} cells, {adata_subset.n_vars} genes")

# -------------------------
# Create output dir
# -------------------------
os.makedirs("figures", exist_ok=True)

# -------------------------
# GLOBAL UMAP (subset only)
# -------------------------
fig, ax = plt.subplots(figsize=(6, 5))

sc.pl.umap(
    adata_subset,
    color="sample",
    size=20,
    ax=ax,
    show=False
)

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_umap.png", dpi=300)
plt.close()

# -------------------------
# PER-SAMPLE UMAP
# -------------------------
samples = adata_subset.obs["sample"].unique()

fig, axes = plt.subplots(1, len(samples), figsize=(5 * len(samples), 5))

# handle single sample case
if len(samples) == 1:
    axes = [axes]

for i, s in enumerate(samples):
    sc.pl.umap(
        adata_subset[adata_subset.obs["sample"] == s],
        color="sample",
        size=20,
        title=f"{s}",
        ax=axes[i],
        show=False
    )

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_perSample_umap.png", dpi=300)
plt.close()

# -------------------------
# Save subset
# -------------------------
adata_subset.write_h5ad(args.output)
print(f"Saved: {args.output}")

# -------------------------
# Summary
# -------------------------
print("\nSummary:")
print(f"  Cells: {adata_subset.n_obs}")
print(f"  Genes: {adata_subset.n_vars}")
print(f"  Samples: {adata_subset.obs['sample'].value_counts().to_dict()}")
