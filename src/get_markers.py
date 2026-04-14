import scanpy as sc
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    adata = sc.read(args.input)

    if "leiden" not in adata.obs:
        raise ValueError("ERROR: 'leiden' not found in adata.obs")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    result = adata.uns["rank_genes_groups"]
    clusters = result["names"].dtype.names

    top_n = 10

    data = []
    for cl in clusters:
        genes = result["names"][cl][:top_n]
        row = [cl] + list(genes)
        data.append(row)

    columns = ["cluster"] + [f"gene_{i+1}" for i in range(top_n)]
    df = pd.DataFrame(data, columns=columns)

    out_csv = f"{args.prefix}_leiden_top_genes.csv"
    df.to_csv(out_csv, index=False)

    print(f"Saved: {out_csv}")

    adata.write(args.output, compression="gzip")

if __name__ == "__main__":
    main()
