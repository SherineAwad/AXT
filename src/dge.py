import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--groupBy", required=True)
parser.add_argument("--N", required=True, type=int)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

celltypes = adata.obs["celltype"].unique()
all_rows = []

for ct in celltypes:
    adata_ct = adata[adata.obs["celltype"] == ct].copy()
    if len(adata_ct.obs["sample"].unique()) < 2:
        print(f"WARNING: {ct} has only one sample, skipping")
        continue
    
    sc.tl.rank_genes_groups(adata_ct, groupby="sample", method="wilcoxon", layer='log1p', use_raw=False)
    
    r = adata_ct.uns["rank_genes_groups"]
    groups = r["names"].dtype.names
    
    for g in groups:
        for i in range(len(r["names"][g])):
            all_rows.append({
                "group": ct,
                "gene": r["names"][g][i],
                "logfoldchanges": r["logfoldchanges"][g][i],
                "pvals_adj": r["pvals_adj"][g][i]
            })

df = pd.DataFrame(all_rows)
df = df[df["pvals_adj"] < 0.05]
df.to_csv(f"{args.prefix}_{args.groupBy}_dge.csv", index=False)

# Get top N genes per celltype sorted by absolute logFC then pval_adj
top_genes_dict = {}
for ct in celltypes:
    ct_df = df[df["group"] == ct]
    ct_df = ct_df.iloc[ct_df["logfoldchanges"].abs().argsort()[::-1]]
    ct_df = ct_df.sort_values(["logfoldchanges", "pvals_adj"], ascending=[False, True])
    top_genes_dict[ct] = ct_df.head(args.N)["gene"].tolist()

ordered_genes = []
for ct in celltypes:
    for gene in top_genes_dict[ct]:
        if gene not in ordered_genes:
            ordered_genes.append(gene)

matrix = []
for gene in ordered_genes:
    row = []
    for ct in celltypes:
        val = df[(df["group"] == ct) & (df["gene"] == gene)]["logfoldchanges"].values
        row.append(val[0] if len(val) > 0 else 0)
    matrix.append(row)

heatmap_df = pd.DataFrame(matrix, index=ordered_genes, columns=celltypes)

gene_max_ct = {}
for gene in ordered_genes:
    gene_data = heatmap_df.loc[gene]
    max_ct = gene_data.abs().idxmax()
    gene_max_ct[gene] = max_ct

sorted_genes = []
for ct in celltypes:
    for gene in ordered_genes:
        if gene_max_ct.get(gene) == ct and gene not in sorted_genes:
            sorted_genes.append(gene)

heatmap_df = heatmap_df.loc[sorted_genes]

plt.figure(figsize=(max(10, len(heatmap_df.columns)), max(8, len(heatmap_df.index) * 0.3)))
sns.heatmap(heatmap_df, cmap="RdBu_r", center=0, annot=False)
plt.title(f"Reg vs nonReg: Log Fold Change per Cell Type")
plt.ylabel("Genes")
plt.xlabel("Cell Type")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_{args.groupBy}_dge_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

adata.write_h5ad(args.output)
