#!/usr/bin/env python3

import argparse
import anndata as ad
import decoupler as dc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--celltype', default='celltype')
parser.add_argument('--prefix', required=True)
parser.add_argument('--organism', default='mouse')
parser.add_argument('--progeny_model', required=True)
parser.add_argument('--top_n', type=int, default=50)
parser.add_argument('--csv', help='Optional CSV with interactions to plot separately')
args = parser.parse_args()

print("Loading...")
adata = ad.read_h5ad(args.input)

if 'liana_res' not in adata.uns:
    raise ValueError("liana_res not found in adata.uns")

full_liana = adata.uns['liana_res'].copy()

liana_df = (
    adata.uns['liana_res']
    .sort_values('lrscore', ascending=False)
    .head(args.top_n)
    .copy()
)

print("Running PROGENY...")

progeny_net = pd.read_csv(args.progeny_model)

progeny_net = progeny_net.rename(
    columns={
        'pathway': 'source',
        'gene': 'target',
        'weight': 'weight'
    }
)

dc.mt.mlm(
    adata,
    progeny_net,
    tmin=1
)

progeny_scores = adata.obsm['score_mlm']

pathway_scores = pd.DataFrame(
    progeny_scores,
    index=adata.obs_names
)

pathway_scores['celltype'] = adata.obs[args.celltype].values

pathway_by_ct = (
    pathway_scores
    .groupby('celltype')
    .mean()
)

def get_progeny_scores_for_interactions(interaction_df):
    results = []
    
    for _, row in interaction_df.iterrows():
        receiver = row['target']
        
        if receiver not in pathway_by_ct.index:
            continue
        
        rec_scores = pathway_by_ct.loc[receiver]
        
        out_row = {
            'sender': row['source'],
            'receiver': receiver,
            'ligand': row['ligand_complex'],
            'receptor': row['receptor_complex'],
            'lrscore': row['lrscore']
        }
        
        for pathway in pathway_by_ct.columns:
            out_row[pathway] = rec_scores[pathway]
        
        results.append(out_row)
    
    return pd.DataFrame(results)

df = get_progeny_scores_for_interactions(liana_df)

os.makedirs("figures", exist_ok=True)

df.to_csv(
    f"{args.prefix}_lr_progeny.csv",
    index=False
)

print(f"Saved {len(df)} interactions")

def make_dotplot(data, filename, title_suffix=""):
    if len(data) == 0:
        print(f"No data for {filename}")
        return
    
    pathway_cols = [
        c for c in data.columns
        if c not in ['sender', 'receiver', 'ligand', 'receptor', 'lrscore']
    ]
    
    plot_df = data.copy()
    
    plot_df['interaction'] = (
        plot_df['sender']
        + ':'
        + plot_df['ligand']
        + ' → '
        + plot_df['receiver']
        + ':'
        + plot_df['receptor']
    )
    
    melted = plot_df.melt(
        id_vars=['interaction'],
        value_vars=pathway_cols,
        var_name='pathway',
        value_name='score'
    )
    
    melted = melted.dropna()
    
    # Create matrix for heatmap
    pivot = melted.pivot(index='interaction', columns='pathway', values='score')
    
    # Plot
    fig, ax = plt.subplots(figsize=(max(10, len(pathway_cols) * 0.6), max(6, len(plot_df) * 0.3)))
    
    # Custom colormap: blue-white-red
    cmap = LinearSegmentedColormap.from_list('progeny', ['blue', 'white', 'red'])
    
    # Get max abs value for symmetric color scale
    max_abs = max(abs(pivot.min().min()), abs(pivot.max().max()))
    if max_abs == 0:
        max_abs = 1
    
    # Create dotplot
    for i, interaction in enumerate(pivot.index):
        for j, pathway in enumerate(pivot.columns):
            score = pivot.loc[interaction, pathway]
            color = cmap((score / max_abs + 1) / 2)  # Normalize to 0-1 for colormap
            
            ax.scatter(j, i, s=100, c=[color], edgecolors='black', linewidth=0.5)
    
    ax.set_xticks(range(len(pathway_cols)))
    ax.set_xticklabels(pathway_cols, rotation=45, ha='right')
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    
    ax.set_xlabel("PROGENY pathway")
    ax.set_ylabel("LIANA interaction")
    ax.set_title(f"LIANA interactions vs PROGENY pathways {title_suffix} ({args.prefix})")
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-max_abs, max_abs))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('PROGENY score')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

if len(df) > 0:
    make_dotplot(df, f"figures/{args.prefix}_lr_progeny_dotplot.png")

if args.csv:
    csv_df = pd.read_csv(args.csv)
    print(f"\nLoaded CSV with {len(csv_df)} interactions")
    
    csv_df['key'] = csv_df['source'] + '|' + csv_df['target'] + '|' + csv_df['ligand_complex'] + '|' + csv_df['receptor_complex']
    full_liana['key'] = full_liana['source'] + '|' + full_liana['target'] + '|' + full_liana['ligand_complex'] + '|' + full_liana['receptor_complex']
    
    csv_keys = set(csv_df['key'])
    csv_liana_rows = full_liana[full_liana['key'].isin(csv_keys)].copy()
    
    print(f"Found {len(csv_liana_rows)} CSV interactions in full LIANA results")
    
    if len(csv_liana_rows) > 0:
        csv_with_scores = get_progeny_scores_for_interactions(csv_liana_rows)
        
        if len(csv_with_scores) > 0:
            make_dotplot(csv_with_scores, f"figures/{args.prefix}_lr_progeny_interest_plot.png", "")
            print(f"Created additional plot for {len(csv_with_scores)} CSV interactions")

adata.write_h5ad(f"{args.prefix}.h5ad")

print("\nDone")
