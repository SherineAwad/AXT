import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, required=True)
parser.add_argument("--pvalue", type=float, default=0.05)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

if "celltype" not in adata.obs.columns:
    raise ValueError("celltype column not found in adata.obs")

sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",
    method="wilcoxon",
    layer="log1p",
    use_raw=False,
    n_genes=None
)

df = sc.get.rank_genes_groups_df(adata, group=None)
df = df[["names", "logfoldchanges", "pvals_adj", "group"]].rename(
    columns={"names": "gene", "group": "celltype"}
)

df.to_csv(f"{args.prefix}_DGEperCelltypeAll.csv", index=False)
df_sig = df[df["pvals_adj"] < args.pvalue].copy()
df_sig.to_csv(f"{args.prefix}_DGEperCelltype_{args.pvalue}.csv", index=False)

if df_sig.empty:
    raise ValueError("No significant genes found")

# Get top N up and down genes per celltype
top_genes = []
for ct in df_sig["celltype"].unique():
    ct_df = df_sig[df_sig["celltype"] == ct]
    
    up = ct_df[ct_df["logfoldchanges"] > 0].sort_values(
        "logfoldchanges", ascending=False
    ).head(args.N)
    
    down = ct_df[ct_df["logfoldchanges"] < 0].sort_values(
        "logfoldchanges", ascending=True
    ).head(args.N)
    
    top_genes.extend(up["gene"].tolist())
    top_genes.extend(down["gene"].tolist())

top_genes = list(set(top_genes))

# Create matrix with ALL logFC values for these genes across ALL celltypes
heatmap_data = []
for gene in top_genes:
    for ct in df_sig["celltype"].unique():
        val = df_sig[(df_sig["gene"] == gene) & (df_sig["celltype"] == ct)]
        if not val.empty:
            heatmap_data.append({
                "gene": gene,
                "celltype": ct,
                "logfoldchanges": val["logfoldchanges"].values[0]
            })

heatmap_df = pd.DataFrame(heatmap_data).pivot(
    index="gene",
    columns="celltype",
    values="logfoldchanges"
)

# Sort genes by which celltype they have strongest effect
celltype_order = heatmap_df.columns.tolist()
gene_order = []
for gene in heatmap_df.index:
    # Find which celltype has highest absolute logFC for this gene
    max_ct = heatmap_df.loc[gene].abs().idxmax()
    gene_order.append((gene, max_ct))

heatmap_df = heatmap_df.loc[[g for g, _ in sorted(gene_order, key=lambda x: celltype_order.index(x[1]))]]

# Plot
plt.figure(figsize=(max(10, len(heatmap_df.columns) * 0.8), max(8, len(heatmap_df) * 0.3)))

sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    cbar_kws={"label": "logFC"},
    square=False
)

plt.title(f"Top {args.N} up/down genes per cell type\n(logFC values from DGE comparing each celltype vs all others)")
plt.ylabel("Genes")
plt.xlabel("Cell Type")
plt.xticks(rotation=45, ha='right')

plt.tight_layout()

os.makedirs("figures", exist_ok=True)
plt.savefig(f"figures/{args.prefix}_celltype_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

adata.write_h5ad(args.output)

print("Done.")
