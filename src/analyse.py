import scanpy as sc
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    adata = sc.read(args.input)

    # Unique names
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Normalise & log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Scale, PCA, neighbours, UMAP
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Figure quality
    sc.set_figure_params(dpi=120, dpi_save=300)

    # Global UMAP colored by sample
    sc.pl.umap(
        adata,
        color='sample',
        size=20,
        save=f"_{args.prefix}.png"
    )

    # Per-sample UMAPs - SAME COORDINATES, just filtered
    for s in adata.obs['sample'].unique():
        sc.pl.umap(
            adata[adata.obs['sample'] == s],
            color='sample',
            size=20,
            title=f"Sample: {s}",
            save=f"_{args.prefix}_{s}.png"
        )

    # Save processed object
    adata.write(args.output)

if __name__ == "__main__":
    main()


