import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def load_sheet(xls, sheet):
    df = pd.read_excel(xls, sheet_name=sheet)

    required = ["name", "p_value", "intersection_size", "query_size"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{sheet} missing columns: {missing}")

    return df


def clean(df, direction):
    df = df[["name", "p_value", "intersection_size", "query_size"]].copy()

    df["p_value"] = df["p_value"].replace(0, np.nextafter(0, 1))
    df["p_value"] = df["p_value"].clip(lower=np.nextafter(0, 1))

    df["neglogp"] = -np.log10(df["p_value"])
    df["gene_ratio"] = df["intersection_size"] / df["query_size"]

    df = df.sort_values("neglogp", ascending=False)
    df = df.drop_duplicates("name")

    df["direction"] = direction
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xls", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    # -------------------------
    # LOAD
    # -------------------------
    up = clean(load_sheet(args.xls, "Up"), "UP")
    down = clean(load_sheet(args.xls, "Down"), "DOWN")

    df = pd.concat([up, down], ignore_index=True)

    if df.empty:
        raise ValueError("No data to plot")

    # -------------------------
    # ORDER - most significant on TOP
    # -------------------------
    df = df.sort_values("neglogp", ascending=False)

    terms = df["name"].drop_duplicates().tolist()
    y_map = {t: i for i, t in enumerate(terms)}
    df["y"] = df["name"].map(y_map)

    # -------------------------
    # DOT SIZE
    # -------------------------
    df["dot_size"] = df["neglogp"] * 15

    # -------------------------
    # PLOT
    # -------------------------
    colors = {"UP": "red", "DOWN": "blue"}

    plt.figure(figsize=(10, max(6, len(terms) * 0.3)))

    for direction in ["UP", "DOWN"]:
        sub = df[df["direction"] == direction]

        plt.scatter(
            sub["gene_ratio"],
            sub["y"],
            s=sub["dot_size"],
            alpha=0.75,
            color=colors[direction],
            edgecolors="black",
            linewidth=0.3
        )

    # -------------------------
    # AXES - invert Y so most significant (smallest y value) appears on TOP
    # -------------------------
    plt.yticks(range(len(terms)), terms, fontsize=8)
    plt.gca().invert_yaxis()  # This puts index 0 (most significant) at the top
    plt.xlabel("Gene Ratio (intersection_size / query_size)")
    plt.ylabel("GO Terms")
    plt.title("GO Enrichment Dotplot (UP vs DOWN)")

    # -------------------------
    # LEGEND 1: COLOR (VERTICAL - NORMAL)
    # -------------------------
    legend_color = [
        Line2D([0], [0], marker='o', color='w', label='UP',
               markerfacecolor='red', markersize=7, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='DOWN',
               markerfacecolor='blue', markersize=7, markeredgecolor='black'),
    ]

    leg1 = plt.legend(
        handles=legend_color,
        title="Direction",
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        ncol=1,
        frameon=True
    )

    # -------------------------
    # LEGEND 2: SIZE (VERTICAL - NORMAL)
    # -------------------------
    size_levels = [5, 10, 20]

    legend_size = [
        Line2D([0], [0],
               marker='o',
               color='w',
               label=f"-log10(p) = {v}",
               markerfacecolor='gray',
               markeredgecolor='black',
               markersize=np.sqrt(v * 15))
        for v in size_levels
    ]

    leg2 = plt.legend(
        handles=legend_size,
        title="Significance (dot size)",
        loc="lower left",
        bbox_to_anchor=(1.02, 0.5),
        ncol=1,
        frameon=True
    )

    plt.gca().add_artist(leg1)

    # -------------------------
    # FINAL
    # -------------------------
    plt.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    out = f"figures/{args.prefix}_GO_enrichment_dotplot.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
