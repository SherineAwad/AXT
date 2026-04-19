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

df.to_csv(f"{args.prefix}_perCelltypeAll.csv", index=False)
df_sig = df[df["pvals_adj"] < args.pvalue].copy()
df_sig.to_csv(f"{args.prefix}_perCelltype_{args.pvalue}.csv", index=False)

if df_sig.empty:
    raise ValueError("No significant genes found")

top_list = []
for ct in df_sig["celltype"].unique():
    ct_df = df_sig[df_sig["celltype"] == ct]
    
    up = ct_df[ct_df["logfoldchanges"] > 0].sort_values(
        "logfoldchanges", ascending=False
    ).head(args.N)
    
    down = ct_df[ct_df["logfoldchanges"] < 0].sort_values(
        "logfoldchanges", ascending=True
    ).head(args.N)
    
    top_list.append(up)
    top_list.append(down)

plot_df = pd.concat(top_list).drop_duplicates(["celltype", "gene"])

heatmap_df = plot_df.pivot(
    index="gene",
    columns="celltype",
    values="logfoldchanges"
).fillna(0)

celltype_order = heatmap_df.columns.tolist()
diag_sort = []
for gene in heatmap_df.index:
    max_ct = heatmap_df.loc[gene].abs().idxmax()
    diag_sort.append((gene, max_ct))

heatmap_df = heatmap_df.loc[[g for g, _ in sorted(diag_sort, key=lambda x: celltype_order.index(x[1]))]]

plt.figure(figsize=(max(8, len(heatmap_df.columns) * 0.5), max(6, len(heatmap_df) * 0.3)))

sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    cbar_kws={"label": "logFC"}
)

plt.title(f"Top {args.N} up and down genes per cell type")
plt.ylabel("Genes")
plt.xlabel("Cell Type")
plt.xticks(rotation=45, ha='right')

plt.tight_layout()

os.makedirs("figures", exist_ok=True)
plt.savefig(f"figures/{args.prefix}_celltype_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

adata.write_h5ad(args.output)

print("Done.")
