import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from gprofiler import GProfiler
import os

parser = argparse.ArgumentParser()
parser.add_argument("--csv", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--organism", default="mmusculus")
parser.add_argument("--top_n", type=int, default=10)
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--source", default="KEGG,REAC")
args = parser.parse_args()

# Source mapping
source_map = {
    "GO": ["GO:BP", "GO:CC", "GO:MF"],
    "BP": ["GO:BP"],
    "CC": ["GO:CC"],
    "MF": ["GO:MF"],
    "KEGG": ["KEGG"],
    "REAC": ["REAC"]
}

source_list = []
for s in args.source.split(","):
    s = s.strip()
    if s in source_map:
        source_list.extend(source_map[s])
    else:
        source_list.append(s)

df = pd.read_csv(args.csv)

if "celltype" not in df.columns:
    df["celltype"] = "all"

if "gene" not in df.columns:
    raise ValueError("CSV must contain a 'gene' column")

gp = GProfiler(return_dataframe=True)

celltypes = df["celltype"].dropna().unique().tolist()

all_results = []

for ct in celltypes:
    print(f"Processing {ct}...")

    ct_df = df[df["celltype"] == ct]

    up_genes = ct_df[ct_df["logfoldchanges"] > 0]["gene"].dropna().tolist()

    if len(up_genes) >= 3:
        try:
            enr_up = gp.profile(
                query=up_genes,
                organism=args.organism,
                sources=source_list,
                user_threshold=args.pvalue
            )

            if not enr_up.empty:
                enr_up["celltype"] = ct
                enr_up["direction"] = "UP"
                enr_up["gene_count"] = enr_up["intersection_size"]
                enr_up["neg_log10_pval"] = -np.log10(enr_up["p_value"])
                all_results.append(enr_up)

        except Exception as e:
            print(f"UP error: {e}")

    down_genes = ct_df[ct_df["logfoldchanges"] < 0]["gene"].dropna().tolist()

    if len(down_genes) >= 3:
        try:
            enr_down = gp.profile(
                query=down_genes,
                organism=args.organism,
                sources=source_list,
                user_threshold=args.pvalue
            )

            if not enr_down.empty:
                enr_down["celltype"] = ct
                enr_down["direction"] = "DOWN"
                enr_down["gene_count"] = enr_down["intersection_size"]
                enr_down["neg_log10_pval"] = -np.log10(enr_down["p_value"])
                all_results.append(enr_down)

        except Exception as e:
            print(f"DOWN error: {e}")

if len(all_results) == 0:
    raise ValueError("No enrichment results found")

combined = pd.concat(all_results, ignore_index=True)
combined.to_csv(f"{args.prefix}_enrichment_{args.source}.csv", index=False)

# =========================================================
# ONLY CHANGE: GSEA BARPLOT WITH SOURCE LABEL
# =========================================================

fig, (ax_up, ax_down) = plt.subplots(1, 2, figsize=(28, 10))

def prepare(direction):
    df_dir = combined[combined["direction"] == direction].copy()

    if df_dir.empty:
        return pd.DataFrame()

    if "source" not in df_dir.columns:
        df_dir["source"] = "UNKNOWN"

    df_grouped = df_dir.groupby(["name", "source"]).agg(
        p_value=("p_value", "min"),
        gene_count=("gene_count", "max")
    ).reset_index()

    df_grouped["score"] = -np.log10(df_grouped["p_value"])
    df_grouped["label"] = df_grouped["name"] + " (" + df_grouped["source"] + ")"

    df_grouped = df_grouped.sort_values("score", ascending=True)

    return df_grouped

def plot(ax, df_plot, title, color):
    if df_plot.empty:
        ax.set_title(f"{title} (no results)")
        return

    top_df = df_plot.tail(args.top_n)

    ax.barh(
        top_df["label"],
        top_df["score"],
        color=color,
        edgecolor="black",
        alpha=0.85
    )

    ax.set_title(title, fontsize=18)
    ax.set_xlabel("-log10(p-value)", fontsize=14)
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='x', labelsize=12)

up_df = prepare("UP")
down_df = prepare("DOWN")

plot(ax_up, up_df, "UP regulated", "firebrick")
plot(ax_down, down_df, "DOWN regulated", "royalblue")

plt.tight_layout()

os.makedirs("figures", exist_ok=True)
plt.savefig(
    f"figures/{args.prefix}_enrichment.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Done.")
