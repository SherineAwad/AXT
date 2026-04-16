import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, required=True)
parser.add_argument("--pvalue", type=float, default=0.05)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# ----------------------------
# GLOBAL DGE
# ----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="sample",
    method="wilcoxon",
    layer="log1p",
    use_raw=False
)

r = adata.uns["rank_genes_groups"]
groups = r["names"].dtype.names

all_rows = []

for g in groups:
    for i in range(len(r["names"][g])):
        all_rows.append({
            "gene": r["names"][g][i],
            "logfoldchanges": r["logfoldchanges"][g][i],
            "pvals_adj": r["pvals_adj"][g][i]
        })

df = pd.DataFrame(all_rows)

df.to_csv(f"{args.prefix}_global_dge.csv", index=False)

# ----------------------------
# TOP GENES
# ----------------------------
up = df[
    (df["logfoldchanges"] > 0) &
    (df["pvals_adj"] < args.pvalue)
].sort_values("logfoldchanges", ascending=False).head(args.N)

down = df[
    (df["logfoldchanges"] < 0) &
    (df["pvals_adj"] < args.pvalue)
].sort_values("logfoldchanges", ascending=True).head(args.N)

plot_df = pd.concat([up, down]).drop_duplicates("gene")

if len(plot_df) == 0:
    raise ValueError("No significant genes found")

# ----------------------------
# HEATMAP = DIRECT FROM CSV (NO RECOMPUTATION)
# ----------------------------
heatmap_df = plot_df.set_index("gene")[["logfoldchanges"]]

# sort for readability
heatmap_df = heatmap_df.sort_values("logfoldchanges", ascending=False)

# ----------------------------
# PLOT
# ----------------------------
plt.figure(figsize=(6, max(6, len(heatmap_df) * 0.25)))

sns.heatmap(
    heatmap_df,
    cmap="RdBu_r",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC (Reg vs non-Reg)"}
)

plt.title("Global DE (Reg vs non-Reg)")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(
    f"figures/{args.prefix}_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

adata.write_h5ad(args.output)
