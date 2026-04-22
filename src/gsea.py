import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from gseapy import prerank
import os

parser = argparse.ArgumentParser()
parser.add_argument("--csv", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--organism", default="Mouse")
parser.add_argument("--top_n", type=int, default=10)
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--source", default="KEGG_2019_Mouse,Reactome_2022_Mouse")
args = parser.parse_args()

df = pd.read_csv(args.csv)
df = df.drop_duplicates(subset=['gene'])

# Sort by absolute logFC, then by p-value
df['abs_logfc'] = df['logfoldchanges'].abs()
df = df.sort_values(['abs_logfc', 'pvals_adj'], ascending=[False, True])

# Create sequential rank
df['rank'] = range(1, len(df) + 1)

# Multiply rank by sign of logFC
df['ranking_score'] = df['rank'] * np.sign(df['logfoldchanges'])

ranked_genes = df.set_index('gene')['ranking_score']

print(f"Ranked {len(ranked_genes)} genes")

gsea_results = prerank(
    rnk=ranked_genes,
    gene_sets=args.source,
    organism=args.organism,
    outdir=None,
    min_size=5,
    max_size=500,
    permutation_num=1000,
    seed=42
)

results = []
for term in gsea_results.results:
    if gsea_results.results[term]['fdr'] < args.pvalue:
        hits = gsea_results.results[term]['hits']
        hits_str = ','.join([str(h) for h in hits[:10]])
        results.append({
            'term': term,
            'nes': gsea_results.results[term]['nes'],
            'fdr': gsea_results.results[term]['fdr'],
            'pval': gsea_results.results[term]['pval'],
            'direction': 'UP' if gsea_results.results[term]['nes'] > 0 else 'DOWN',
            'genes': hits_str
        })

df_results = pd.DataFrame(results)

if df_results.empty:
    raise ValueError(f"No significant results at FDR < {args.pvalue}")

df_results.to_csv(f"{args.prefix}_gsea.csv", index=False)

df_up = df_results[df_results['direction'] == 'UP'].head(args.top_n)
df_down = df_results[df_results['direction'] == 'DOWN'].head(args.top_n)
df_plot = pd.concat([df_up, df_down])

fig, ax = plt.subplots(figsize=(10, max(6, len(df_plot) * 0.4)))

colors = ['#D62728'] * len(df_up) + ['#1F77B4'] * len(df_down)

ax.barh(df_plot['term'], df_plot['nes'], color=colors, alpha=0.8)
ax.axvline(0, color='black', linewidth=0.5)
ax.set_xlabel('NES')
ax.set_ylabel('Pathway')
ax.set_title('GSEA: UP vs DOWN')

plt.tight_layout()
os.makedirs("figures", exist_ok=True)
plt.savefig(f"figures/{args.prefix}_gsea.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Found {len(df_results)} pathways | UP: {len(df_up)} | DOWN: {len(df_down)}")
