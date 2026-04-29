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
parser.add_argument("--input", required=True)
parser.add_argument("--genes", required=True)
parser.add_argument("--prefix", default="gene_expression")
parser.add_argument("--topn", type=int, default=20)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
adata = sc.read_h5ad(args.input)

if "rank_genes_groups" not in adata.uns:
    raise ValueError("No DGE found")

# ----------------------------
# LOAD SECRETED GENES
# ----------------------------
with open(args.genes) as f:
    secreted = set(g.strip() for g in f if g.strip())

secreted = secreted.intersection(set(adata.var_names))

# ----------------------------
# BUILD DGE TABLE
# ----------------------------
dfs = []
for g in adata.uns["rank_genes_groups"]["names"].dtype.names:
    df_g = sc.get.rank_genes_groups_df(adata, group=g)
    df_g["group"] = g
    dfs.append(df_g)

df = pd.concat(dfs, ignore_index=True)
df["neglog10_padj"] = -np.log10(df["pvals_adj"] + 1e-300)

# ----------------------------
# TOP N SELECTION
# ----------------------------
df_unique = df.drop_duplicates("names")

df_top = (
    df_unique.sort_values(
        by=["pvals_adj", "logfoldchanges"],
        ascending=[True, False]
    )
    .head(args.topn)
)

top_genes = set(df_top["names"])

# ----------------------------
# SECRETED GENES
# ----------------------------
secreted_in_dge = set(df[df["names"].isin(secreted)]["names"])

# ----------------------------
# FINAL GENES
# ----------------------------
genes_plot = sorted(top_genes.union(secreted_in_dge))

print(f"Top genes: {len(top_genes)}")
print(f"Secreted genes: {len(secreted_in_dge)}")
print(f"Total plotted genes: {len(genes_plot)}")

df_plot = df[df["names"].isin(genes_plot)]

# ============================================================
# VOLCANO PLOT (FIXED COLOR LOGIC)
# ============================================================
plt.figure(figsize=(14, 10))

colors = []
for _, r in df.iterrows():
    fc = r["logfoldchanges"]

    if fc > 1:
        colors.append("cornflowerblue")
    elif fc < -1:
        colors.append("lightcoral")
    else:
        colors.append("lightgrey")

plt.scatter(
    df["logfoldchanges"],
    df["neglog10_padj"],
    s=5,
    alpha=0.15,
    c=colors
)

highlight_genes = set(secreted)
df_highlight = df[df["names"].isin(highlight_genes)]

plt.scatter(
    df_highlight["logfoldchanges"],
    df_highlight["neglog10_padj"],
    s=30,
    color="black",
    edgecolors="black"
)

for _, r in df_highlight.iterrows():
    plt.annotate(
        r["names"],
        xy=(r["logfoldchanges"], r["neglog10_padj"]),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=8,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            lw=0.6,
            shrinkA=0,
            shrinkB=0
        )
    )

plt.axhline(-np.log10(0.05), ls="--", color="blue", alpha=0.5)
plt.axvline(0, color="black", alpha=0.3)

plt.xlabel("log Fold Change")
plt.ylabel("-log10 adj p-value")
plt.title(f"Volcano: Top {args.topn} + Secreted Genes")

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_volcano.png", dpi=300)
plt.close()

# ============================================================
# FEATURE PLOTS (UNCHANGED)
# ============================================================
if "X_umap" in adata.obsm:
    for gene in genes_plot:
        if gene not in adata.var_names:
            continue

        try:
            sc.pl.umap(
                adata,
                color=gene,
                layer="log1p",
                show=False,
                title=f"{gene}",
                color_map="viridis",
                splitby="sample"
            )

            plt.savefig(
                f"figures/{args.prefix}_feature_{gene}.png",
                dpi=300,
                bbox_inches="tight"
            )
            plt.close()

        except:
            pass

# ============================================================
# DOTPLOT (UNCHANGED)
# ============================================================
if len(genes_plot) > 0:
    try:
        sc.pl.dotplot(
            adata,
            var_names=genes_plot,
            groupby="sample",
            layer="log1p",
            show=False,
            dendrogram=False,
            swap_axes=True,
            standard_scale="var",
            figsize=(max(8, len(genes_plot) * 0.4), 6)
        )

        plt.xticks(rotation=90, fontsize=9)
        plt.yticks(fontsize=10)

        plt.tight_layout()

        plt.savefig(
            f"figures/{args.prefix}_dotplot.png",
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()

    except:
        pass

# ============================================================
# VIOLIN (UNCHANGED)
# ============================================================
for gene in genes_plot:
    if gene not in adata.var_names:
        continue

    try:
        df_expr = []
        samples = list(adata.obs["sample"].unique())

        for s in samples:
            vals = adata[adata.obs["sample"] == s, gene].layers["log1p"]
            vals = vals.toarray().flatten() if hasattr(vals, "toarray") else vals.flatten()

            for v in vals:
                df_expr.append({"sample": s, "expr": v})

        df_expr = pd.DataFrame(df_expr)

        sample_order = samples[::-1]

        palette = {
            samples[0]: "lightcoral",
            samples[1] if len(samples) > 1 else samples[0]: "cornflowerblue"
        }

        plt.figure()

        ax = sns.violinplot(
            data=df_expr,
            x="sample",
            y="expr",
            order=sample_order,
            palette=palette,
            inner=None
        )

        for artist in ax.collections:
            artist.set_edgecolor("black")

        plt.xticks(rotation=45)
        plt.tight_layout()

        plt.savefig(
            f"figures/{args.prefix}_violin_{gene}.png",
            dpi=300
        )
        plt.close()

    except:
        pass

print("DONE -> figures/")
