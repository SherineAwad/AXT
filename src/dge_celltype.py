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
parser.add_argument("--reference", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

all_results = []

# ----------------------------
# PER CELLTYPE DE
# ----------------------------
for ct in adata.obs["celltype"].unique():

    adata_ct = adata[adata.obs["celltype"] == ct].copy()

    if adata_ct.n_obs < 10:
        continue

    sc.tl.rank_genes_groups(
        adata_ct,
        groupby="sample",
        reference=args.reference,
        method="wilcoxon",
        layer="log1p",
        use_raw=False
    )

    df_ct = sc.get.rank_genes_groups_df(adata_ct, group=None)

    df_ct = df_ct[["names", "logfoldchanges", "pvals_adj"]].rename(
        columns={"names": "gene"}
    )

    df_ct["celltype"] = ct

    all_results.append(df_ct)

df = pd.concat(all_results, ignore_index=True)

df[df["pvals_adj"] < args.pvalue].to_csv(
    f"{args.prefix}_perCelltype_dge.csv", index=False)

# ----------------------------
# FILTER
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()

sig = sig.sort_values("logfoldchanges", key=lambda x: x.abs(), ascending=False)
sig = sig.drop_duplicates(["celltype", "gene"])

# ----------------------------
# TOP N
# ----------------------------
top_list = []

for ct in sig["celltype"].unique():
    ct_df = sig[sig["celltype"] == ct]

    up = ct_df[ct_df["logfoldchanges"] > 0].sort_values(
        "logfoldchanges", ascending=False
    ).head(args.N)

    down = ct_df[ct_df["logfoldchanges"] < 0].sort_values(
        "logfoldchanges", ascending=True
    ).head(args.N)

    top_list.append(up)
    top_list.append(down)

plot_df = pd.concat(top_list).drop_duplicates(["celltype", "gene"])

if plot_df.empty:
    raise ValueError("No significant genes found")

# ----------------------------
# HEATMAP
# ----------------------------
heatmap_df = plot_df.pivot(
    index="gene",
    columns="celltype",
    values="logfoldchanges"
).fillna(0)

# keep previous improved ordering (unchanged)
gene_order = heatmap_df.mean(axis=1).sort_values(ascending=False).index
heatmap_df = heatmap_df.loc[gene_order]

# ----------------------------
# PLOT
# ----------------------------
plt.figure(figsize=(max(6, len(heatmap_df.columns)), max(6, len(heatmap_df) * 0.3)))

sns.heatmap(
    heatmap_df,
    cmap="RdBu_r",
    center=0,
    cbar_kws={"label": "logFC (Reg vs nonReg)"}
)

plt.title("Per-celltype DE: Reg vs nonReg")
plt.ylabel("Genes")
plt.xlabel("celltype")

plt.tight_layout()

os.makedirs("figures", exist_ok=True)
plt.savefig(f"figures/{args.prefix}_celltype_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

adata.write_h5ad(args.output)
