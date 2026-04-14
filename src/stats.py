#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

os.makedirs("figures", exist_ok=True)

# Create counts table (celltype x sample)
ct_counts = pd.crosstab(adata.obs["celltype"], adata.obs["sample"])
ct_counts["Total"] = ct_counts.sum(axis=1)
ct_counts.loc["Total"] = ct_counts.sum(axis=0)

# Save counts CSV (NO FRACTIONS)
ct_counts.to_csv(f"{args.prefix}_cell_counts.csv")

# For plot: calculate fractions per sample
ratios_df = ct_counts.drop("Total").drop(columns="Total").div(ct_counts.drop("Total").drop(columns="Total").sum(axis=0), axis=1)

# Stacked bar plot
ax = ratios_df.T.plot(kind='bar', stacked=True, figsize=(12, 6))
plt.ylabel("Fraction of cells")
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample")
plt.xticks(rotation=45, ha='right')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_cell_ratios.png", dpi=600, bbox_inches='tight')
plt.close()

print(f"Saved: {args.prefix}_cell_counts.csv")
print(f"Saved: figures/{args.prefix}_cell_ratios.png")
