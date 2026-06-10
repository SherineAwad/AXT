#!/usr/bin/env python3

import argparse
import anndata as ad
import decoupler as dc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--celltype', default='celltype')
parser.add_argument('--prefix', required=True)
parser.add_argument('--organism', default='mouse')
parser.add_argument('--progeny_model', required=True)
parser.add_argument('--top_n', type=int, default=50)
args = parser.parse_args()

print("Loading...")
adata = ad.read_h5ad(args.input)

if 'liana_res' not in adata.uns:
    raise ValueError("liana_res not found in adata.uns")

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

results = []

for _, row in liana_df.iterrows():

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

df = pd.DataFrame(results)

os.makedirs("figures", exist_ok=True)

df.to_csv(
    f"{args.prefix}_lr_progeny.csv",
    index=False
)

print(f"Saved {len(df)} interactions")

if len(df) > 0:

    pathway_cols = [
        c for c in df.columns
        if c not in ['sender', 'receiver', 'ligand', 'receptor', 'lrscore']
    ]

    plot_df = df.copy()

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

    scaler = MinMaxScaler()

    melted['size'] = scaler.fit_transform(
        np.abs(melted[['score']])
    ) * 300 + 50

    plt.figure(
        figsize=(
            max(8, len(pathway_cols) * 0.8),
            max(6, len(plot_df) * 0.4)
        )
    )

    plt.scatter(
        x=pd.Categorical(
            melted['pathway'],
            categories=pathway_cols
        ).codes,
        y=pd.Categorical(
            melted['interaction'],
            categories=plot_df['interaction'].unique()
        ).codes,
        s=melted['size'],
        c=melted['score'],
        cmap='coolwarm',
        edgecolors='black'
    )

    plt.xticks(
        range(len(pathway_cols)),
        pathway_cols,
        rotation=45,
        ha='right'
    )

    plt.yticks(
        range(len(plot_df['interaction'].unique())),
        plot_df['interaction'].unique(),
        fontsize=8
    )

    plt.xlabel("PROGENY pathway")
    plt.ylabel("LIANA interaction")
    plt.title(
        f"LIANA interactions vs PROGENY pathways ({args.prefix})"
    )

    plt.colorbar(label="PROGENY score")

    plt.tight_layout()

    plt.savefig(
        f"figures/{args.prefix}_lr_progeny_dotplot.png",
        dpi=300,
        bbox_inches='tight'
    )

    plt.close()

adata.write_h5ad(args.output)

print("Done")
