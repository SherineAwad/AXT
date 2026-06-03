#!/usr/bin/env python3

import scanpy as sc
import cellrank as cr
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # -------------------------
    # args
    # -------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cluster_key", default="celltype")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    # -------------------------
    # load data
    # -------------------------
    adata = sc.read_h5ad(args.input)
    print(f"[DEBUG] Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # -------------------------
    # preprocessing
    # -------------------------
    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"]
        print("[DEBUG] Using log1p layer")

    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata)
        print("[DEBUG] PCA computed")

    sc.pp.neighbors(adata)

    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    # -------------------------
    # PAGA - Compute BEFORE CellRank
    # -------------------------
    print("[DEBUG] Computing PAGA...")
    sc.tl.paga(adata, groups=args.cluster_key)

    # -------------------------
    # CellRank (NO root, NO pseudotime)
    # -------------------------
    print("[DEBUG] Building ConnectivityKernel...")

    k = cr.kernels.ConnectivityKernel(adata)
    k.compute_transition_matrix()

    print("[DEBUG] Initializing GPCCA...")
    gpi = cr.estimators.GPCCA(k)

    print("[DEBUG] Schur decomposition...")
    gpi.compute_schur()

    if args.cluster_key not in adata.obs:
        raise ValueError(f"{args.cluster_key} not found in adata.obs")

    print("[DEBUG] Computing macrostates...")
    gpi.compute_macrostates(cluster_key=args.cluster_key)

    # push macrostates to adata
    if hasattr(gpi, "macrostates"):
        adata.obs["macrostates"] = gpi.macrostates

    print("[DEBUG] Predicting terminal states...")
    gpi.predict_terminal_states()

    # -------------------------
    # REQUIRED ORDER FIX
    # -------------------------
    print("[DEBUG] Computing fate probabilities...")
    gpi.compute_fate_probabilities()

    print("[DEBUG] Computing lineage drivers...")
    gpi.compute_lineage_drivers()

    # -------------------------
    # PLOTS
    # -------------------------
    print("[DEBUG] UMAP macrostates...")

    if "macrostates" in adata.obs and adata.obs["macrostates"].nunique() > 1:
        sc.pl.umap(
            adata,
            color="macrostates",
            save=f"_{args.prefix}_macrostates.png",
            show=False
        )
    else:
        print("\n" + "="*80)
        print("BIOLOGY INTERPRETATION WARNING (NOT A CODE ERROR)")
        print("="*80)
        print(
            "Only ONE macrostate detected.\n"
            "This indicates no separable metastable structure under Connectivity-only CellRank.\n\n"
            "This is a BIOLOGICAL / DATA STRUCTURE limitation, NOT a software error.\n"
        )
        print("="*80 + "\n")

    # Only plot fate probabilities if multiple lineages exist
    if hasattr(gpi, "fate_probabilities") and gpi.fate_probabilities.shape[1] > 1:
        print("[DEBUG] Fate probabilities plot...")
        cr.pl.aggregate_fate_probabilities(
            adata,
            basis="umap",
            cluster_key=args.cluster_key,
            save=f"_{args.prefix}_fates.png",
            show=False
        )
    else:
        print("[DEBUG] Skipping fate probabilities plot - only 1 lineage detected")

    print("[DEBUG] PAGA plot...")
    sc.pl.paga(
        adata,
        color=args.cluster_key,
        save=f"_{args.prefix}_paga.png",
        show=False
    )

    # -------------------------
    # SUMMARY
    # -------------------------
    print("\n" + "="*50)
    print("TERMINAL STATES")
    print("="*50)

    terminal_states = getattr(gpi, "terminal_states", None)

    if terminal_states is not None:
        for s in list(terminal_states):
            cells = adata.obs[adata.obs["macrostates"] == s]
            if len(cells) > 0:
                dom = cells[args.cluster_key].value_counts().index[0]
                print(f"{s}: {dom} ({len(cells)})")

    print("\nDONE")
