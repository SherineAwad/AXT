import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input h5ad file")
parser.add_argument("--genes", required=True, help="Text file with genes (one per line)")
parser.add_argument("--prefix", default="gene_expression", help="Prefix for output files")
args = parser.parse_args()

# ----------------------------
# CREATE OUTPUT DIR
# ----------------------------
os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
print(f"Loading data from {args.input}")
adata = sc.read_h5ad(args.input)

# ----------------------------
# CHECK FOR LOG1P LAYER
# ----------------------------
if "log1p" not in adata.layers:
    raise ValueError("log1p layer not found in adata")

# ----------------------------
# READ GENES
# ----------------------------
print(f"Reading genes from {args.genes}")
with open(args.genes, 'r') as f:
    genes = [line.strip() for line in f if line.strip()]

# ----------------------------
# VALIDATE GENES
# ----------------------------
missing_genes = [g for g in genes if g not in adata.var_names]
if missing_genes:
    print(f"Warning: {len(missing_genes)} genes not found: {missing_genes}")
    genes = [g for g in genes if g in adata.var_names]
    
if not genes:
    raise ValueError("No valid genes found")

print(f"Plotting {len(genes)} genes")

# ----------------------------
# 1. DOTPLOT
# ----------------------------
print("Creating dotplot...")
fig, ax = plt.subplots(figsize=(12, max(3, len(genes) * 0.4)))

sc.pl.dotplot(
    adata,
    var_names=genes,
    groupby="sample",
    layer="log1p",
    use_raw=False,
    ax=ax,
    show=False,
    standard_scale="var",
    cmap="coolwarm",
    dot_max=0.5,
    title=f"Gene Expression: {', '.join(genes[:5])}{'...' if len(genes)>5 else ''}"
)

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_dotplot.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------
# 2. VIOLIN PLOTS (one per gene)
# ----------------------------
print("Creating violin plots...")
for gene in genes:
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Extract data for this gene
    gene_data = []
    for sample in adata.obs["sample"].unique():
        sample_mask = adata.obs["sample"] == sample
        expr = adata[sample_mask, gene].layers["log1p"]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        
        for val in expr:
            gene_data.append({"sample": sample, "expression": val})
    
    df_gene = pd.DataFrame(gene_data)
    
    # Custom colors
    samples = df_gene["sample"].unique()
    colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD", "#98D8C8", "#F7B05E"]
    palette = {sample: colors[i % len(colors)] for i, sample in enumerate(samples)}
    
    sns.violinplot(data=df_gene, x="sample", y="expression", ax=ax, 
                  palette=palette, inner="box", linewidth=1.5)
    
    ax.set_title(f"{gene}", fontsize=14, fontweight="bold")
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("log1p Expression", fontsize=12)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"figures/{args.prefix}_violin_{gene}.png", dpi=300, bbox_inches="tight")
    plt.close()

# ----------------------------
# 3. BOXPLOT WITH JITTER (one per gene)
# ----------------------------
print("Creating boxplot with jitter...")
for gene in genes:
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Extract data for this gene
    gene_data = []
    for sample in adata.obs["sample"].unique():
        sample_mask = adata.obs["sample"] == sample
        expr = adata[sample_mask, gene].layers["log1p"]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        
        for val in expr:
            gene_data.append({"sample": sample, "expression": val})
    
    df_gene = pd.DataFrame(gene_data)
    
    # Custom colors
    samples = df_gene["sample"].unique()
    colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD", "#98D8C8", "#F7B05E"]
    palette = {sample: colors[i % len(colors)] for i, sample in enumerate(samples)}
    
    # Boxplot with light grey
    sns.boxplot(data=df_gene, x="sample", y="expression", ax=ax, 
               color="#E8E8E8", width=0.6, linewidth=1.5, boxprops=dict(edgecolor="#555555"))
    
    # Stripplot with sample colors
    for sample in samples:
        sample_data = df_gene[df_gene["sample"] == sample]
        ax.scatter(sample_data["sample"], sample_data["expression"], 
                  alpha=0.6, s=20, color=palette[sample], label=sample)
    
    ax.set_title(f"{gene}", fontsize=14, fontweight="bold")
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("log1p Expression", fontsize=12)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"figures/{args.prefix}_boxplot_jitter_{gene}.png", dpi=300, bbox_inches="tight")
    plt.close()

# ----------------------------
# SUMMARY
# ----------------------------
print("\n" + "="*50)
print("SUMMARY")
print("="*50)
print(f"Input file: {args.input}")
print(f"Number of cells: {adata.n_obs}")
print(f"Genes plotted: {len(genes)}")
print(f"Samples found: {adata.obs['sample'].unique().tolist()}")
print(f"Figures saved to: figures/")
print(f"  - {args.prefix}_dotplot.png")
print(f"  - {args.prefix}_violin_GENE.png ({len(genes)} files)")
print(f"  - {args.prefix}_boxplot_jitter_GENE.png ({len(genes)} files)")
print("="*50)
