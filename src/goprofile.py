import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from gprofiler import GProfiler
import os

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--csv", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--organism", default="mmusculus")
parser.add_argument("--top_n", type=int, default=100)
parser.add_argument("--plot_n", type=int, default=10)
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--source", default="KEGG,REAC")
args = parser.parse_args()

# ----------------------------
# SOURCE MAPPING
# ----------------------------
source_map = {
    "GO": ["GO:BP", "GO:CC", "GO:MF"],
    "BP": ["GO:BP"],
    "CC": ["GO:CC"],
    "MF": ["GO:MF"],
    "KEGG": ["KEGG"],
    "REAC": ["REAC"],
    "CORUM": ["CORUM"],
    "WP": ["WP"]
}

source_list = []
for s in args.source.split(","):
    s = s.strip()
    if s in source_map:
        source_list.extend(source_map[s])
    else:
        source_list.append(s)

# ----------------------------
# LOAD DATA
# ----------------------------
df = pd.read_csv(args.csv)

# enforce required columns
required = ["gene", "logfoldchanges", "pvals_adj"]
for col in required:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

if "celltype" not in df.columns:
    df["celltype"] = "all"

gp = GProfiler(return_dataframe=True)

celltypes = df["celltype"].dropna().unique().tolist()
all_results = []

# ----------------------------
# MAIN LOOP
# ----------------------------
for ct in celltypes:
    print(f"Processing {ct}...")

    ct_df = df[df["celltype"] == ct].copy()

    # ============================
    # UP REGULATED
    # ============================
    up_df = ct_df[ct_df["logfoldchanges"] > 0].copy()
    up_df = up_df.sort_values("pvals_adj", ascending=True)

    up_genes = up_df["gene"].dropna().head(args.top_n).tolist()

    if len(up_genes) >= 5:
        try:
            enr_up = gp.profile(
                query=up_genes,
                organism=args.organism,
                sources=source_list
            )

            if not enr_up.empty:
                enr_up["celltype"] = ct
                enr_up["direction"] = "UP"
                enr_up["gene_count"] = enr_up["intersection_size"]
                enr_up["neg_log10_pval"] = -np.log10(
                    enr_up["p_value"].clip(lower=1e-300)
                )
                all_results.append(enr_up)

        except Exception as e:
            print(f"UP error: {e}")

    # ============================
    # DOWN REGULATED
    # ============================
    down_df = ct_df[ct_df["logfoldchanges"] < 0].copy()
    down_df = down_df.sort_values("pvals_adj", ascending=True)

    down_genes = down_df["gene"].dropna().head(args.top_n).tolist()

    if len(down_genes) >= 5:
        try:
            enr_down = gp.profile(
                query=down_genes,
                organism=args.organism,
                sources=source_list
            )

            if not enr_down.empty:
                enr_down["celltype"] = ct
                enr_down["direction"] = "DOWN"
                enr_down["gene_count"] = enr_down["intersection_size"]
                enr_down["neg_log10_pval"] = -np.log10(
                    enr_down["p_value"].clip(lower=1e-300)
                )
                all_results.append(enr_down)

        except Exception as e:
            print(f"DOWN error: {e}")

# ----------------------------
# CHECK RESULTS
# ----------------------------
if len(all_results) == 0:
    raise ValueError("No enrichment results found")

combined = pd.concat(all_results, ignore_index=True)

# filter AFTER enrichment
combined = combined[combined["p_value"] <= args.pvalue]

# save table
os.makedirs("results", exist_ok=True)
combined.to_csv(
    f"results/{args.prefix}_enrichment_{args.source}.csv",
    index=False
)

# ----------------------------
# PLOTTING
# ----------------------------
fig, (ax_up, ax_down) = plt.subplots(1, 2, figsize=(20, 10))

def prepare(df, direction):
    df_dir = df[df["direction"] == direction].copy()

    if df_dir.empty:
        return pd.DataFrame()

    df_dir = df_dir.sort_values("neg_log10_pval", ascending=True)

    df_dir["label"] = df_dir["name"] + " (" + df_dir["source"] + ")"

    return df_dir.tail(args.plot_n)

def plot(ax, df_plot, title, color):
    if df_plot.empty:
        ax.set_title(f"{title} (no results)")
        return

    ax.barh(
        df_plot["label"],
        df_plot["neg_log10_pval"],
        color=color,
        edgecolor="black",
        alpha=0.85
    )

    ax.set_title(title, fontsize=16)
    ax.set_xlabel("-log10(p-value)", fontsize=12)
    ax.tick_params(axis='y', labelsize=9)

up_plot = prepare(combined, "UP")
down_plot = prepare(combined, "DOWN")

plot(ax_up, up_plot, "UP regulated pathways", "firebrick")
plot(ax_down, down_plot, "DOWN regulated pathways", "royalblue")

plt.tight_layout()

os.makedirs("figures", exist_ok=True)
plt.savefig(
    f"figures/{args.prefix}_{args.source}_enrichment.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Done.")
