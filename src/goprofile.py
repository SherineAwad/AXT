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

# Parse source
source_list = []
for s in args.source.split(","):
    s = s.strip()
    if s in source_map:
        source_list.extend(source_map[s])
    else:
        source_list.append(s)

# Read CSV
df = pd.read_csv(args.csv)

# Initialize g:Profiler
gp = GProfiler(return_dataframe=True)

# Get unique cell types
celltypes = df["celltype"].unique()

# Run enrichment
all_results = []

for ct in celltypes:
    print(f"Processing {ct}...")
    
    up_genes = df[(df["celltype"] == ct) & (df["logfoldchanges"] > 0)]["gene"].tolist()
    
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
            print(f"  Error: {e}")
    
    down_genes = df[(df["celltype"] == ct) & (df["logfoldchanges"] < 0)]["gene"].tolist()
    
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
            print(f"  Error: {e}")

if not all_results:
    raise ValueError("No enrichment results found")

# Save CSV
combined = pd.concat(all_results, ignore_index=True)
combined.to_csv(f"{args.prefix}_enrichment_{args.source}.csv", index=False)

# Create figure with 3 columns: UP | DOWN | LEGEND
fig = plt.figure(figsize=(42, 16))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.3])

ax_up = fig.add_subplot(gs[0])
ax_down = fig.add_subplot(gs[1])
ax_legend = fig.add_subplot(gs[2])

def get_marker(gene_count):
    if gene_count <= 5:
        return 'o'      # circle
    elif gene_count <= 10:
        return 's'      # square
    elif gene_count <= 15:
        return '^'      # triangle
    else:
        return 'D'      # diamond

# UP panel - UNCHANGED
dir_df = combined[combined["direction"] == "UP"]
if not dir_df.empty:
    top_terms = []
    for ct in celltypes:
        ct_df = dir_df[dir_df["celltype"] == ct]
        if not ct_df.empty:
            top = ct_df.sort_values("p_value").head(args.top_n)
            top_terms.append(top)
    top_terms_df = pd.concat(top_terms)
    
    for ct in top_terms_df["celltype"].unique():
        ct_data = top_terms_df[top_terms_df["celltype"] == ct]
        for _, row in ct_data.iterrows():
            marker = get_marker(row["gene_count"])
            ax_up.scatter(ct, row["name"], 
                         s=25, marker=marker, c=[row["neg_log10_pval"]], 
                         cmap="Reds", vmin=0, vmax=5, alpha=0.8,
                         edgecolors='black', linewidth=0.5)
    
    ax_up.set_title("UP regulated", fontsize=18)
    ax_up.set_xlabel("Cell Type", fontsize=16)
    ax_up.set_ylabel("Pathway", fontsize=16)
    ax_up.tick_params(axis='x', rotation=45, labelsize=14)
    ax_up.tick_params(axis='y', labelsize=12)

# DOWN panel - UNCHANGED
dir_df = combined[combined["direction"] == "DOWN"]
if not dir_df.empty:
    top_terms = []
    for ct in celltypes:
        ct_df = dir_df[dir_df["celltype"] == ct]
        if not ct_df.empty:
            top = ct_df.sort_values("p_value").head(args.top_n)
            top_terms.append(top)
    top_terms_df = pd.concat(top_terms)
    
    for ct in top_terms_df["celltype"].unique():
        ct_data = top_terms_df[top_terms_df["celltype"] == ct]
        for _, row in ct_data.iterrows():
            marker = get_marker(row["gene_count"])
            ax_down.scatter(ct, row["name"], 
                           s=25, marker=marker, c=[row["neg_log10_pval"]], 
                           cmap="Reds", vmin=0, vmax=5, alpha=0.8,
                           edgecolors='black', linewidth=0.5)
    
    ax_down.set_title("DOWN regulated", fontsize=18)
    ax_down.set_xlabel("Cell Type", fontsize=16)
    ax_down.set_ylabel("Pathway", fontsize=16)
    ax_down.tick_params(axis='x', rotation=45, labelsize=14)
    ax_down.tick_params(axis='y', labelsize=12)

# ========================
# LEGEND PANEL - FIXED
# ========================
ax_legend.axis('off')

# 1. Colorbar at top (taller)
cax = ax_legend.inset_axes([0.2, 0.55, 0.6, 0.4])  # x, y, width, height (taller)
sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0, 5))
sm.set_array([])
cbar = plt.colorbar(sm, cax=cax, orientation='vertical')
cbar.set_label('-log10(p-value)', fontsize=12)

# 2. Tight box for gene counts at bottom
box_x, box_y, box_width, box_height = 0.15, 0.05, 0.7, 0.4
rect = plt.Rectangle((box_x, box_y), box_width, box_height, 
                     fill=True, facecolor='white', edgecolor='black', linewidth=1.5)
ax_legend.add_patch(rect)

# Title inside box
ax_legend.text(box_x + box_width/2, box_y + box_height - 0.05, 'Gene count', 
               fontsize=11, ha='center', fontweight='bold')

# Markers and labels (tight vertical stack inside box)
marker_info = [
    ('o', '1-5 genes'),
    ('s', '6-10 genes'),
    ('^', '11-15 genes'),
    ('D', '16+ genes')
]

y_start = box_y + box_height - 0.12
for i, (marker, label) in enumerate(marker_info):
    ax_legend.scatter(box_x + 0.15, y_start - i*0.07, s=50, marker=marker, 
                      c='gray', edgecolors='black')
    ax_legend.text(box_x + 0.35, y_start - i*0.07, label, fontsize=10, va='center')

plt.tight_layout()
os.makedirs("figures", exist_ok=True)
plt.savefig(f"figures/{args.prefix}_enrichment_{args.source}_dotplot.png", dpi=300, bbox_inches="tight")
plt.close()

print("Done.")
