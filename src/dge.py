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
parser.add_argument("--reference", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# ----------------------------
# GLOBAL DE
# ----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="sample",
    reference=args.reference,
    method="wilcoxon",
    layer="log1p",
    use_raw=False
)

# ----------------------------
# CLEAN EXTRACTION (NO r["..."][0])
# ----------------------------
df = sc.get.rank_genes_groups_df(adata, group=None)

df = df[["names", "logfoldchanges", "pvals_adj"]].rename(
    columns={"names": "gene"}
)

# ----------------------------
# SAVE CSV
# ----------------------------

df[df["pvals_adj"] < args.pvalue].to_csv(
    f"{args.prefix}_global_dge.csv", index=False)
# ----------------------------
# FILTER
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()

# remove duplicates (keep strongest signal)
sig = sig.sort_values("logfoldchanges", key=lambda x: x.abs(), ascending=False)
sig = sig.drop_duplicates("gene")

# ----------------------------
# TOP N UP / DOWN
# ----------------------------
up = sig[sig["logfoldchanges"] > 0].sort_values(
    "logfoldchanges", ascending=False
).head(args.N)

down = sig[sig["logfoldchanges"] < 0].sort_values(
    "logfoldchanges", ascending=True
).head(args.N)

plot_df = pd.concat([up, down]).drop_duplicates("gene")

if plot_df.empty:
    raise ValueError("No significant genes found")

# ----------------------------
# HEATMAP
# ----------------------------
heatmap_df = plot_df.set_index("gene")[["logfoldchanges"]]

heatmap_df = heatmap_df.reindex(
    heatmap_df["logfoldchanges"].abs().sort_values(ascending=False).index
)

# ----------------------------
# PLOT
# ----------------------------
plt.figure(figsize=(6, max(6, len(heatmap_df) * 0.3)))

sns.heatmap(
    heatmap_df,
    cmap="RdBu_r",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC"}
)

plt.title(f"Global DE: Reg vs non Reg")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

adata.write_h5ad(args.output)
