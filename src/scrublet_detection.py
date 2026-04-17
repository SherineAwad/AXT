import scrublet as scr
import scanpy as sc
import numpy as np
from scipy import sparse
import argparse
import os
import matplotlib.pyplot as plt

# -------------------------
# Args
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--prefix', required=True)
parser.add_argument('--threshold', type=float, default=0.5)

args = parser.parse_args()

# -------------------------
# Read AnnData object
# -------------------------
adata = sc.read_h5ad(args.input)

# Fix observation names warning
if not adata.obs_names.is_unique:
    adata.obs_names_make_unique()

# -------------------------
# Get expression matrix (dense for Scrublet)
# -------------------------
print("Preparing data for Scrublet...")
if sparse.issparse(adata.X):
    X = adata.X.toarray()
else:
    X = adata.X

# -------------------------
# Run Scrublet doublet detection
# -------------------------
print(f"Running Scrublet on {adata.n_obs} cells...")
scrub = scr.Scrublet(X, expected_doublet_rate=0.06)  # 6% expected for 31k cells

# Simulate doublets and get scores
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

# -------------------------
# Apply threshold
# -------------------------
# Scrublet uses its own threshold, but we can apply voter_thresh equivalent
threshold = args.threshold
predicted_doublets_final = doublet_scores > threshold

# -------------------------
# Figures
# -------------------------
os.makedirs("figures", exist_ok=True)

# Histogram of doublet scores
fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(doublet_scores, bins=50, edgecolor='black')
ax.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold = {threshold}')
ax.set_xlabel('Doublet Score')
ax.set_ylabel('Frequency')
ax.set_title(f'{args.prefix} - Scrublet Doublet Scores')
ax.legend()
fig.savefig(f"figures/{args.prefix}_scrublet_scores.png", dpi=300)
plt.close(fig)

# -------------------------
# Save results
# -------------------------
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets_final

print(f"Detected {np.sum(predicted_doublets_final)} doublets ({np.mean(predicted_doublets_final)*100:.1f}%)")

adata.write(args.output, compression="gzip")

print("Scrublet doublet detection completed and results saved.")
