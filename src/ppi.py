#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as ad
import omnipath as op
import decoupler as dc
import networkx as nx
import textwrap

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--prefix', required=True)
parser.add_argument('--organism', default='mouse')
parser.add_argument('--top_n', type=int, default=50)
parser.add_argument('--path_len', type=int, default=None,
                    help='Only consider TFs with path length <= this value (e.g., 1 for direct interactions)')
args = parser.parse_args()

# Load data
print("Loading data...")
adata = ad.read_h5ad(args.input)
if 'liana_res' not in adata.uns:
    raise ValueError("liana_res not found")
liana_df = adata.uns['liana_res'].sort_values('lrscore', ascending=False).head(args.top_n).copy()

# Build undirected graph (case‑insensitive)
print("Loading OmniPath...")
ppi = op.interactions.OmniPath().get(genesymbols=True)
G = nx.Graph()
for _, r in ppi.iterrows():
    s, t = r['source_genesymbol'], r['target_genesymbol']
    if pd.notna(s) and pd.notna(t):
        G.add_edge(s.upper(), t.upper())
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Load TFs (CollecTRI)
print("Loading TFs...")
tf_net = dc.op.collectri(organism=args.organism, remove_complexes=False, license='academic')
TFs = {tf.upper() for tf in tf_net['source']}
G.add_nodes_from(TFs)
print(f"Total TFs: {len(TFs)}")

def first_gene(complex_str):
    if pd.isna(complex_str):
        return None
    for sep in ['_', ':', '-']:
        if sep in complex_str:
            return complex_str.split(sep)[0].strip()
    return complex_str.strip()

results = []
for _, row in liana_df.iterrows():
    sender_ct = row['source']
    receiver_ct = row['target']
    ligand = row['ligand_complex']
    receptor_complex = row['receptor_complex']
    lrscore = row['lrscore']

    rec_gene = first_gene(receptor_complex)
    if rec_gene is None:
        results.append({
            'sender_ct': sender_ct, 'receiver_ct': receiver_ct, 'ligand': ligand,
            'receptor': receptor_complex, 'lrscore': lrscore,
            'best_tf': '', 'best_len': np.nan, 'integration_score': lrscore
        })
        continue

    rec_upper = rec_gene.upper()
    if rec_upper not in G:
        results.append({
            'sender_ct': sender_ct, 'receiver_ct': receiver_ct, 'ligand': ligand,
            'receptor': receptor_complex, 'lrscore': lrscore,
            'best_tf': '', 'best_len': np.nan, 'integration_score': lrscore
        })
        continue

    # Find shortest path to any TF within the path_len limit
    best_len = np.inf
    best_tf = ''
    for tf in TFs:
        if tf not in G:
            continue
        try:
            path = nx.shortest_path(G, rec_upper, tf)
            length = len(path) - 1
            if args.path_len is not None and length > args.path_len:
                continue
            if length < best_len:
                best_len = length
                best_tf = tf
            elif length == best_len and best_tf > tf:  # tie‑break alphabetically
                best_tf = tf
        except nx.NetworkXNoPath:
            continue

    if best_tf == '':
        # No TF found within limit
        results.append({
            'sender_ct': sender_ct, 'receiver_ct': receiver_ct, 'ligand': ligand,
            'receptor': receptor_complex, 'lrscore': lrscore,
            'best_tf': '', 'best_len': np.nan, 'integration_score': lrscore
        })
    else:
        integration_score = lrscore + 1 / (best_len + 1)
        results.append({
            'sender_ct': sender_ct, 'receiver_ct': receiver_ct, 'ligand': ligand,
            'receptor': receptor_complex, 'lrscore': lrscore,
            'TF': best_tf, 'Path Length': best_len, 'integration_score': integration_score
        })

df = pd.DataFrame(results)

# Save outputs
os.makedirs("figures", exist_ok=True)
csv_path = f"{args.prefix}_liana_omnipath_tf_best.csv"
df.to_csv(csv_path, index=False)
adata.uns['liana_omnipath_tf'] = df
adata.write_h5ad(args.output)
print(f"Saved H5AD → {args.output}")
print(f"Saved CSV → {csv_path}")

# Plot only pairs that have a best_tf
print("Plotting...")
df_plot = df[df['best_tf'] != ''].sort_values("integration_score", ascending=True)
if len(df_plot) == 0:
    print("No ligand–receptor pairs found with TF reachable within the path length limit.")
else:
    labels = []
    for _, row in df_plot.iterrows():
        lr_label = f"{row['sender_ct']}.{row['ligand']} → {row['receiver_ct']}.{row['receptor']}"
        label = f"{lr_label}\nBest TF: {row['best_tf']} (len {int(row['best_len'])})"
        wrapped = "\n".join(textwrap.wrap(label, width=70))
        labels.append(wrapped)
    height = max(0.5, len(df_plot) * 0.4)
    plt.figure(figsize=(12, height))
    plt.barh(labels, df_plot['integration_score'], color='tab:green')
    plt.xlabel("Integration score = lrscore + 1/(path_len+1)")
    plt.title(f"Ligand-Receptor → best TF ({args.prefix})")
    plt.tight_layout()
    fig_path = f"figures/{args.prefix}_liana_omnipath_tf_best.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved plot → {fig_path}")
print("Done.")
