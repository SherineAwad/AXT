#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import liana as li
import plotnine as p9


def make_key(df):
    return (
        df["source"].astype(str)
        + "->"
        + df["target"].astype(str)
        + "|"
        + df["ligand_complex"].astype(str)
        + "-"
        + df["receptor_complex"].astype(str)
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample1", required=True)
    parser.add_argument("--sample2", required=True)
    parser.add_argument("--prefix", required=True)

    # MODE SELECTION
    parser.add_argument("--ligand", default=None, help="Filter by SOURCE cell type")
    parser.add_argument("--receptor", default=None, help="Filter by TARGET cell type")

    parser.add_argument("--lrscore", type=float, default=0.90)
    parser.add_argument("--specificity_rank", type=float, default=0.01)

    parser.add_argument("--common", choices=["true", "false"], default="false")

    args = parser.parse_args()

    if args.ligand and args.receptor:
        raise ValueError("Use ONLY --ligand OR --receptor, not both.")

    if not args.ligand and not args.receptor:
        raise ValueError("You must provide either --ligand or --receptor.")

    os.makedirs("figures", exist_ok=True)

    # LOAD
    s1 = pd.read_csv(args.sample1, sep=None, engine="python")
    s2 = pd.read_csv(args.sample2, sep=None, engine="python")

    # FILTER MODE
    if args.ligand:
        s1 = s1[s1["source"].astype(str).str.lower() == args.ligand.lower()].copy()
        s2 = s2[s2["source"].astype(str).str.lower() == args.ligand.lower()].copy()
        mode = "ligand"

    if args.receptor:
        s1 = s1[s1["target"].astype(str).str.lower() == args.receptor.lower()].copy()
        s2 = s2[s2["target"].astype(str).str.lower() == args.receptor.lower()].copy()
        mode = "receptor"

    # FILTERS
    s1 = s1[
        (s1["lrscore"] >= args.lrscore)
        & (s1["specificity_rank"] <= args.specificity_rank)
    ].copy()

    s2 = s2[
        (s2["lrscore"] >= args.lrscore)
        & (s2["specificity_rank"] <= args.specificity_rank)
    ].copy()

    if s1.empty:
        print("No interactions in sample1 after filtering.")
        return

    # KEYS
    s1["key"] = make_key(s1)
    s2["key"] = make_key(s2)

    s2_keys = set(s2["key"])

    # COMMON / UNIQUE
    if args.common == "true":
        s1 = s1[s1["key"].isin(s2_keys)].copy()
        suffix = "common"
    else:
        s1 = s1[~s1["key"].isin(s2_keys)].copy()
        suffix = os.path.splitext(os.path.basename(args.sample1))[0]

    if s1.empty:
        print("No interactions after comparison step.")
        return

    # PREP
    s1["specificity"] = s1["specificity_rank"]
    s1 = s1.sort_values("lrscore", ascending=False)

    top_n = min(20, len(s1))

    # -------------------------
    # STRICT CSV EXPORT
    # -------------------------
    out_csv = f"{args.prefix}_{mode}_{suffix}_dotplot.csv"

    cols = [
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "lr_means",
        "cellphone_pvals",
        "expr_prod",
        "scaled_weight",
        "lr_logfc",
        "spec_weight",
        "lrscore",
        "specificity_rank",
        "magnitude_rank",
    ]

    s1.head(top_n).reindex(columns=cols).to_csv(out_csv, index=False)

    print("Saved CSV:", out_csv)

    # DOTPLOT (UNCHANGED)
    plot = li.pl.dotplot(
        liana_res=s1,
        colour="lrscore",
        size="specificity",
        inverse_size=True,
        orderby="lrscore",
        orderby_ascending=False,
        orderby_absolute=True,
        top_n=top_n,
        size_range=(0.5, 4)
    )

    plot = (
        plot
        + p9.theme_bw(base_size=14)
        + p9.scale_color_cmap("RdBu_r")
        + p9.theme(
            axis_text_x=p9.element_text(angle=90),
            figure_size=(12, 7)
        )
    )

    out_file = f"figures/{args.prefix}_{mode}_{suffix}_dotplot.png"
    plot.save(out_file, dpi=300, verbose=False)

    print("Saved:", out_file)
    print("Rows:", len(s1))


if __name__ == "__main__":
    main()
