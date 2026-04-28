import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from scipy.sparse import issparse
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)

LAYER = "log1p"   #  ONLY INPUT LAYER USED

def get_X(adata):
    X = adata.layers[LAYER]
    return X.toarray() if issparse(X) else X

def safe_subset_copy(adata, idx):
    try:
        return adata[idx].copy()
    except Exception:
        X = adata.X[idx]
        X = X.toarray() if issparse(X) else np.array(X)

        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()

        obsm = {}
        if hasattr(adata, "obsm"):
            for k, v in adata.obsm.items():
                try:
                    obsm[k] = v[idx]
                except:
                    obsm[k] = v

        layers = {}
        if hasattr(adata, "layers"):
            for k, v in adata.layers.items():
                try:
                    layers[k] = v[idx] if not issparse(v) else v[idx].toarray()
                except:
                    pass

        return sc.AnnData(X=X, obs=obs, var=var, obsm=obsm, layers=layers)

def run_rf(sub, sample_groups, hvg_n=2000, max_cells=3000, random_state=42):

    # -------------------------
    # split using the two samples
    # -------------------------
    sample1_name, sample2_name = sample_groups[0], sample_groups[1]
    sample1 = sub[sub.obs["sample"] == sample1_name].copy()
    sample2 = sub[sub.obs["sample"] == sample2_name].copy()

    if sample1.n_obs < 10 or sample2.n_obs < 10:
        print("Skipping: insufficient cells")
        return None, None

    # -------------------------
    # subsample
    # -------------------------
    rng = np.random.default_rng(random_state)

    if sample1.n_obs > max_cells:
        idx = rng.choice(sample1.n_obs, max_cells, replace=False)
        sample1 = sample1[idx].copy()

    if sample2.n_obs > max_cells:
        idx = rng.choice(sample2.n_obs, max_cells, replace=False)
        sample2 = sample2[idx].copy()

    # -------------------------
    # gene filtering (LAYER ONLY)
    # -------------------------
    X_s1 = get_X(sample1)
    keep = np.var(X_s1, axis=0) > 0

    sample1 = sample1[:, keep].copy()
    sample2 = sample2[:, keep].copy()

    # HVG on sample1 using log1p layer
    try:
        sc.pp.highly_variable_genes(
            sample1,
            n_top_genes=hvg_n,
            flavor="seurat",
            subset=True,
            layer=LAYER
        )
    except Exception:
        X_tmp = get_X(sample1)
        var_idx = np.argsort(np.var(X_tmp, axis=0))[-hvg_n:]
        sample1 = sample1[:, var_idx].copy()

    sample2 = sample2[:, sample1.var_names].copy()

    # -------------------------
    # matrices (STRICT layer)
    # -------------------------
    X_s1 = np.nan_to_num(get_X(sample1))
    X_s2 = np.nan_to_num(get_X(sample2))

    X_train = np.vstack([X_s1, X_s2])
    y_train = np.hstack([
        np.ones(len(X_s1)),   # first sample = 1
        np.zeros(len(X_s2))   # second sample = 0
    ])

    # -------------------------
    # Random Forest
    # -------------------------
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("rf", RandomForestClassifier(
            n_estimators=200,
            random_state=random_state,
            n_jobs=-1,
            class_weight="balanced"
        ))
    ])

    model.fit(X_train, y_train)

    # -------------------------
    # scoring
    # -------------------------
    score_s1 = model.predict_proba(X_s1)[:, 1]
    score_s2 = model.predict_proba(X_s2)[:, 1]

    sample1.obs["fidelity"] = score_s1
    sample2.obs["fidelity"] = score_s2

    return sample1, sample2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True, help="Prefix for output files (e.g., AXT, control, treatment)")
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--max_cells", type=int, default=3000)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if LAYER not in adata.layers:
        raise ValueError(f"Missing required layer: {LAYER}")

    if "celltype" not in adata.obs:
        raise ValueError("Missing celltype column")

    if "sample" not in adata.obs:
        raise ValueError("Missing sample column")

    # Get the two sample groups dynamically
    sample_groups = sorted(adata.obs["sample"].unique())
    if len(sample_groups) != 2:
        raise ValueError(f"Expected exactly 2 samples, found {len(sample_groups)}")

    print(f"Comparing: {sample_groups[0]} vs {sample_groups[1]}")
    
    adata.obs["fidelity"] = np.nan

    for ct in adata.obs["celltype"].unique():
        print(f"Processing cell type: {ct}")
        
        idx = np.where(adata.obs["celltype"] == ct)[0]
        sub = safe_subset_copy(adata, idx)

        s1, s2 = run_rf(
            sub,
            sample_groups,
            hvg_n=args.hvg_genes,
            max_cells=args.max_cells
        )

        if s1 is not None:
            adata.obs.loc[s1.obs_names, "fidelity"] = s1.obs["fidelity"].values
        if s2 is not None:
            adata.obs.loc[s2.obs_names, "fidelity"] = s2.obs["fidelity"].values

    # -------------------------
    # plot with dynamic filename in figures/ directory
    # -------------------------
    os.makedirs("figures", exist_ok=True)
    
    data = [
        adata.obs.loc[adata.obs["sample"] == sample_groups[0], "fidelity"].dropna().values,
        adata.obs.loc[adata.obs["sample"] == sample_groups[1], "fidelity"].dropna().values
    ]

    # Skip plotting if no data
    if len(data[0]) == 0 or len(data[1]) == 0:
        print(f"Warning: No fidelity scores to plot for {args.prefix}")
    else:
        fig, ax = plt.subplots(figsize=(5,4))

        vp = ax.violinplot(data, showmedians=True)
        for b in vp["bodies"]:
            b.set_alpha(0.7)

        ax.set_xticks([1,2])
        ax.set_xticklabels([sample_groups[0], sample_groups[1]])
        ax.set_ylabel(f"{sample_groups[0]} identity probability")
        ax.set_ylim(0, 1)

        plt.tight_layout()
        figure_filename = os.path.join("figures", f"{args.prefix}_fidelity_violin.png")
        plt.savefig(figure_filename, dpi=300)
        print(f"Saved figure: {figure_filename}")
        plt.close()

    adata.write(args.output)
    print(f"Done: {args.output}")

if __name__ == "__main__":
    main()
