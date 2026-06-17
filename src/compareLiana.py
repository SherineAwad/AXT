#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import liana as li
import plotnine as p9


def make_key(df):
    return (
        df["source"].astype(str) + "->" +
        df["target"].astype(str) + "|" +
        df["ligand_complex"].astype(str) + "-" +
        df["receptor_complex"].astype(str)
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample1", required=True)
    parser.add_argument("--sample2", required=True)
    parser.add_argument("--prefix", required=True)

    parser.add_argument("--lrscore", type=float, default=0.80,
                        help="Minimum LR score threshold")

    parser.add_argument("--specificity_rank", type=float, default=0.05,
                        help="Maximum specificity rank threshold (lower = more specific)")

    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    # -----------------------
    # LOAD DATA
    # -----------------------
    s1 = pd.read_csv(args.sample1, sep=None, engine="python")
    s2 = pd.read_csv(args.sample2, sep=None, engine="python")

    s1["key"] = make_key(s1)
    s2["key"] = make_key(s2)

    # -----------------------
    # FILTER SAMPLE1
    # -----------------------
    s1 = s1[
        (s1["lrscore"] >= args.lrscore) &
        (s1["specificity_rank"] <= args.specificity_rank)
    ].copy()

    # REMOVE anything also in sample2
    s1 = s1[~s1["key"].isin(set(s2["key"]))].copy()

    if s1.empty:
        print("No sample1-specific interactions found.")
        return

    # -----------------------
    # FORMAT FOR LIANA
    # -----------------------
    s1["interaction"] = (
        s1["source"] + "→" +
        s1["target"] + " | " +
        s1["ligand_complex"] + "-" +
        s1["receptor_complex"]
    )

    # LIANA fields
    s1["interaction_stat"] = s1["lrscore"]
    s1["specificity"] = s1["specificity_rank"]

    # sort: strongest first
    s1 = s1.sort_values("interaction_stat", ascending=False)

    top_n = min(10, len(s1))
    s1 = s1.head(top_n)

    # -----------------------
    # LIANA DOTPLOT
    # -----------------------
    plot = li.pl.dotplot(
        liana_res=s1,
        colour='interaction_stat',   # LR strength
        size='specificity',          # specificity rank
        inverse_size=True,           # low rank = big dot
        orderby='interaction_stat',
        orderby_ascending=False,
        orderby_absolute=True,
        top_n=top_n,
        size_range=(0.5, 4)
    )

    # -----------------------
    # STYLE
    # -----------------------
    plot = (
        plot
        + p9.theme_bw(base_size=14)
        + p9.scale_color_cmap('RdBu_r')
        + p9.theme(
            axis_text_x=p9.element_text(angle=90),
            figure_size=(11, 6)
        )
    )

    # -----------------------
    # SAVE
    # -----------------------
    outpath = f"figures/{args.prefix}_liana_lr_specificity_dotplot.png"
    plot.save(outpath, dpi=300)

    print(f"Saved: {outpath}")



if __name__ == "__main__":
    main()
