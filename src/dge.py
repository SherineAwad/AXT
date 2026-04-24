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
parser.add_argument("--N", type=int, required=True)  # heatmap size
parser.add_argument("--K", type=int, required=True)  # volcano labels
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--reference", required=True)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
adata = sc.read_h5ad(args.input)

# ----------------------------
# GLOBAL DGE (UNCHANGED)
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

# Save full DGE (all genes)
full_out = f"{args.prefix}_dgeAll.csv"
df.to_csv(full_out, index=False)

# ----------------------------
# FILTER SIGNIFICANT
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()
sig_out = f"{args.prefix}_dge_{args.pvalue}.csv"
sig.to_csv(sig_out, index=False)

if sig.empty:
    raise ValueError("No significant genes found")

# ============================================================
#  HEATMAP (TOP N BY P-VALUE + DIRECTION)
# ============================================================

top_up_hm = sig[sig["logfoldchanges"] > 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, False]
).head(args.N)

top_down_hm = sig[sig["logfoldchanges"] < 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, True]
).head(args.N)

top_heat = pd.concat([top_up_hm, top_down_hm])

heatmap_df = top_heat.set_index("gene")[["logfoldchanges"]]

plt.figure(figsize=(4, max(3, len(heatmap_df) * 0.12)))


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

# ============================================================
#  VOLCANO (LABEL ONLY TOP K)
# ============================================================

df["neg_log10_pval"] = -np.log10(df["pvals_adj"] + 1e-300)
df["significant"] = df["pvals_adj"] < args.pvalue

# top K for labeling ONLY
top_up_lab = sig[sig["logfoldchanges"] > 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, False]
).head(args.K)

top_down_lab = sig[sig["logfoldchanges"] < 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, True]
).head(args.K)

label_df = pd.concat([top_up_lab, top_down_lab])

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

# highlight top K labeled genes
plt.scatter(
    label_df["logfoldchanges"],
    -np.log10(label_df["pvals_adj"] + 1e-300),
    s=40,
    color="red"
)

# ----------------------------
# LABELS (ONLY TOP K NOW)
# ----------------------------
for _, row in label_df.iterrows():
    x = row["logfoldchanges"]
    y = -np.log10(row["pvals_adj"] + 1e-300)

    x_offset = 2.5 if x > 0 else -2.5
    y_offset = 3.0

    plt.annotate(
        row["gene"],
        xy=(x, y),
        xytext=(x + x_offset, y + y_offset),
        fontsize=5,
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
# SAVE H5AD (UNCHANGED)
# ----------------------------
adata.write_h5ad(args.output)
