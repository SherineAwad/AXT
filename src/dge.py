import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import argparse
import os

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, required=True)
parser.add_argument("--K", type=int, required=True)
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--reference", required=True)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
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

df = sc.get.rank_genes_groups_df(adata, group=None)

df = df[["names", "logfoldchanges", "pvals_adj"]].rename(
    columns={"names": "gene"}
)

# ----------------------------
# SAVE CSVs
# ----------------------------
df[df["pvals_adj"] < args.pvalue].to_csv(
    f"{args.prefix}_global_dge.csv", index=False
)

df.to_csv(
    f"{args.prefix}_DGEperSampleAll.csv", index=False
)

df[df["pvals_adj"] < args.pvalue].to_csv(
    f"{args.prefix}_DGEperSample_{args.pvalue}.csv", index=False
)

# ----------------------------
# FILTER SIGNIFICANT
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()

if sig.empty:
    raise ValueError("No significant genes found")

# ----------------------------
# HEATMAP (TOP N)
# ----------------------------
top_heat = sig.loc[
    sig["logfoldchanges"].abs().sort_values(ascending=False).index
].head(args.N)

heatmap_df = top_heat.set_index("gene")[["logfoldchanges"]]

heatmap_df = heatmap_df.reindex(
    heatmap_df["logfoldchanges"].abs().sort_values(ascending=False).index
)

plt.figure(figsize=(6, max(6, len(heatmap_df) * 0.3)))

sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC"}
)

plt.title(f"Global DE: Reg vs {args.reference}")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------
# VOLCANO PLOT (TOP K UP/DOWN)
# ----------------------------

df["neg_log10_pval"] = -np.log10(df["pvals_adj"] + 1e-300)
df["significant"] = df["pvals_adj"] < args.pvalue

top_up = sig[sig["logfoldchanges"] > 0].nlargest(args.K, "logfoldchanges")
top_down = sig[sig["logfoldchanges"] < 0].nsmallest(args.K, "logfoldchanges")

plot_df = pd.concat([top_up, top_down])

plt.figure(figsize=(8, 6))

# background
plt.scatter(
    df.loc[~df["significant"], "logfoldchanges"],
    df.loc[~df["significant"], "neg_log10_pval"],
    s=10,
    alpha=0.2,
    color="grey"
)

# significant
plt.scatter(
    df.loc[df["significant"], "logfoldchanges"],
    df.loc[df["significant"], "neg_log10_pval"],
    s=10,
    alpha=0.5,
    color="lightgrey"
)

# top UP
plt.scatter(
    top_up["logfoldchanges"],
    -np.log10(top_up["pvals_adj"] + 1e-300),
    s=40,
    color="red"
)

# top DOWN
plt.scatter(
    top_down["logfoldchanges"],
    -np.log10(top_down["pvals_adj"] + 1e-300),
    s=40,
    color="blue"
)

# ----------------------------
# LABELS (ONLY CHANGE: LONGER POINTER OFFSET)
# ----------------------------
for _, row in plot_df.iterrows():
    x = row["logfoldchanges"]
    y = -np.log10(row["pvals_adj"] + 1e-300)

    # LONGER OFFSETS (this is the ONLY change)
    x_offset = 2.5 if x > 0 else -2.5
    y_offset = 3.0

    plt.annotate(
        row["gene"],
        xy=(x, y),
        xytext=(x + x_offset, y + y_offset),
        fontsize=4,
        ha="left" if x > 0 else "right",
        arrowprops=dict(arrowstyle="-", lw=0.5, color="black")
    )

plt.axhline(-np.log10(args.pvalue), linestyle="--")
plt.axvline(0, linestyle="--")

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adjusted p-value")
plt.title(f"Volcano Plot: Reg vs {args.reference}")

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_volcano.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------
# SAVE H5AD
# ----------------------------
adata.write_h5ad(args.output)
