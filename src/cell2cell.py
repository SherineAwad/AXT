import scanpy as sc
import pandas as pd
import liana as li
import argparse
import os
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--database", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--prefix", required=True)

parser.add_argument("--N", type=int, default=10)
parser.add_argument("--logFC", type=float, default=0.5)

args = parser.parse_args()


def run_liana(adata, database, celltype_key):
    li.mt.rank_aggregate(
        adata,
        groupby=celltype_key,
        resource_name=database,
        layer="log1p",
        use_raw=False,
        verbose=True,
        n_perms=1000
    )
    return adata.uns["liana_res"]


print("📥 Loading data...")
adata = sc.read_h5ad(args.input)

res = run_liana(adata, args.database, args.celltype)

# SAVE RAW OUTPUT (UNCHANGED)
res.to_csv(f"{args.prefix}_liana.csv", index=False)

# FILTER (UNCHANGED)
filtered = res[res["lr_logfc"].abs() >= args.logFC].copy()

if filtered.empty:
    raise ValueError("❌ No interactions after filtering")

# TOP N (UNCHANGED)
top_lr = filtered.sort_values("lrscore", ascending=False).head(args.N)

# LABELS (UNCHANGED)
top_lr["ligand_gene"] = top_lr["source"] + "." + top_lr["ligand_complex"]
top_lr["receptor_gene"] = top_lr["target"] + "." + top_lr["receptor_complex"]

plot_df = top_lr.copy()

# ----------------------------
# 🔥 FIX: NUMERIC AXES (NO EMPTY SPACE)
# ----------------------------
ligands = plot_df["ligand_gene"].unique()
receptors = plot_df["receptor_gene"].unique()

lig_map = {l: i for i, l in enumerate(ligands)}
rec_map = {r: i for i, r in enumerate(receptors)}

plot_df["x"] = plot_df["ligand_gene"].map(lig_map)
plot_df["y"] = plot_df["receptor_gene"].map(rec_map)

# ----------------------------
# PLOT (DENSE DOTPLOT)
# ----------------------------
os.makedirs("figures", exist_ok=True)

plt.rcParams.update({"font.size": 12})


plt.figure(figsize=(len(ligands) * 0.6, len(receptors) * 0.5))

plt.scatter(
    plot_df["x"],
    plot_df["y"],
    s=plot_df["lrscore"] * 500,   # scale dots properly
    c=plot_df["lr_logfc"],
    cmap="coolwarm",
    edgecolors="black"
)

# ticks mapped back to labels
plt.xticks(range(len(ligands)), ligands, rotation=45, ha="right")
plt.yticks(range(len(receptors)), receptors)

plt.xlabel("Ligand")
plt.ylabel("Receptor")
plt.title(f"Top ligand–receptor interactions ({args.prefix})")

plt.colorbar(label="lr_logfc")

plt.tight_layout()

plt.savefig(
    f"figures/{args.prefix}_liana_dotplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

adata.write_h5ad(args.output)
