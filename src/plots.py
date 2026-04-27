import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input h5ad file")
parser.add_argument("--genes", required=True, help="Text file with genes (one per line)")
parser.add_argument("--prefix", default="gene_expression", help="Prefix for output files")
parser.add_argument("--minLogFC", type=float, default=0, help="Minimum absolute log fold change to plot (default: 0 = plot all)")
args = parser.parse_args()

# ----------------------------
# VOLCANO PLOT
# ----------------------------
print("Creating volcano plot...")

# Read h5ad (already has DGE results)
adata = sc.read_h5ad(args.input)

# Read genes to highlight
with open(args.genes, 'r') as f:
    highlight_genes = [line.strip() for line in f if line.strip()]

# Extract existing DGE results (NOT re-running)
de_results = sc.get.rank_genes_groups_df(adata, group=None)
de_results = de_results.drop_duplicates(subset='names')

# Calculate -log10 p-value
de_results['neg_log10_pval'] = -np.log10(de_results['pvals_adj'])

# Filter genes below minLogFC cutoff (for plotting background)
de_results_filtered = de_results[abs(de_results['logfoldchanges']) >= args.minLogFC]

print(f"Filtered from {len(de_results)} to {len(de_results_filtered)} genes (|logFC| >= {args.minLogFC})")

# Create volcano plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot filtered background genes
ax.scatter(de_results_filtered['logfoldchanges'],
           de_results_filtered['neg_log10_pval'],
           alpha=0.5, s=10, color='gray')

# ONLY highlight and label genes from --genes file that also pass filter
highlight_df = de_results_filtered[de_results_filtered['names'].isin(highlight_genes)]

for idx, row in highlight_df.iterrows():
    ax.scatter(row['logfoldchanges'],
               row['neg_log10_pval'],
               color='red', s=50)
    ax.annotate(row['names'],
                (row['logfoldchanges'], row['neg_log10_pval']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8)

ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-Log10 Adjusted P-value')
ax.set_title(f'{args.prefix} - Volcano Plot (|logFC| >= {args.minLogFC})')
ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_volcano.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved volcano plot to figures/{args.prefix}_volcano.png")

# ----------------------------
# FEATURE PLOTS (one plot per gene, split by sample)
# ----------------------------
print("Creating feature plots...")

# Get unique samples
samples = adata.obs['sample'].unique()

for gene in highlight_genes:
    if gene in adata.var_names:
        fig, axes = plt.subplots(1, len(samples), figsize=(5*len(samples), 5))
        if len(samples) == 1:
            axes = [axes]
        
        for ax, sample in zip(axes, samples):
            sample_adata = adata[adata.obs['sample'] == sample]
            sc.pl.umap(sample_adata, color=gene, ax=ax, show=False, title=f"{sample}\n{gene}")
        
        plt.suptitle(f"{gene} expression across samples")
        plt.tight_layout()
        plt.savefig(f"figures/{args.prefix}_{gene}_featureplot.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved feature plot for {gene}")

# ----------------------------
# VIOLIN PLOTS for each gene
# ----------------------------
print("Creating violin plots...")

for gene in highlight_genes:
    if gene in adata.var_names:
        sc.pl.violin(adata, gene, groupby='sample', show=False, 
                     save=f"_{args.prefix}_{gene}_violin.png")
        print(f"Saved violin plot for {gene}")

# ----------------------------
# DOTPLOT for all genes
# ----------------------------
print("Creating dotplot...")

sc.pl.dotplot(adata, highlight_genes, groupby='sample', 
              show=False, save=f"_{args.prefix}_dotplot.png")

print(f"Saved dotplot to figures/{args.prefix}_dotplot.png")
