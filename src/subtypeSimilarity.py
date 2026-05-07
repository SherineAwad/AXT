import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("--resolution", type=float, default=0.5, help="Leiden resolution for subclustering")
    parser.add_argument("--min_cells", type=int, default=30, help="Minimum cells to keep a subcluster")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if "celltype" not in adata.obs:
        raise ValueError("Missing celltype column")
    
    if "sample" not in adata.obs:
        raise ValueError("Missing sample column")

    samples = adata.obs["sample"].unique()
    if len(samples) != 2:
        raise ValueError(f"Expected exactly 2 samples, found {len(samples)}")
    
    sample_a, sample_b = samples[0], samples[1]

    print(f"Sample A: {sample_a}")
    print(f"Sample B: {sample_b}")

    celltypes = adata.obs["celltype"].unique()
    
    os.makedirs("figures", exist_ok=True)
    
    all_results = []
    similarity_scores = {}
    
    for ct in celltypes:
        print(f"\nProcessing: {ct}")
        
        # Subset to this cell type
        ct_mask = adata.obs["celltype"] == ct
        sub = adata[ct_mask].copy()
        
        if sub.n_obs < 50:
            print(f"  Skipping: too few cells ({sub.n_obs})")
            similarity_scores[ct] = np.nan
            continue
        
        # Normalize and PCA if not already done
        if "X_pca" not in sub.obsm:
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
            sub = sub[:, sub.var.highly_variable]
            sc.pp.scale(sub)
            sc.tl.pca(sub, n_comps=min(30, sub.n_obs - 1))
        
        # Subcluster with Leiden
        sc.pp.neighbors(sub, n_neighbors=15)
        sc.tl.leiden(sub, resolution=args.resolution, key_added="subcluster")
        
        n_subclusters = len(sub.obs["subcluster"].unique())
        print(f"  Found {n_subclusters} subclusters")
        
        # Count cells per subcluster per sample
        subcluster_data = []
        
        for subc in sub.obs["subcluster"].unique():
            subc_mask = sub.obs["subcluster"] == subc
            subc_cells = sub[subc_mask]
            
            count_a = (subc_cells.obs["sample"] == sample_a).sum()
            count_b = (subc_cells.obs["sample"] == sample_b).sum()
            total = count_a + count_b
            
            if total >= args.min_cells:
                prop_a = count_a / total
                prop_b = count_b / total
                
                subcluster_data.append({
                    "celltype": ct,
                    "subcluster": str(subc),
                    "count_a": count_a,
                    "count_b": count_b,
                    "total": total,
                    "prop_a": prop_a,
                    "prop_b": prop_b,
                    "abs_diff": abs(prop_a - prop_b)
                })
        
        if not subcluster_data:
            print(f"  No subclusters passed min_cells threshold")
            similarity_scores[ct] = np.nan
            continue
        
        df_sub = pd.DataFrame(subcluster_data)
        
        # Calculate similarity score
        avg_diff = df_sub["abs_diff"].mean()
        similarity = 1 - avg_diff
        similarity_scores[ct] = similarity
        
        print(f"  Average difference across subclusters: {avg_diff:.3f}")
        print(f"  Similarity score: {similarity:.3f}")
        
        df_sub = df_sub.sort_values("abs_diff", ascending=False)
        all_results.append(df_sub)
        
        # Plot for this cell type
        fig, ax = plt.subplots(figsize=(10, max(5, len(df_sub) * 0.4)))
        
        y_pos = np.arange(len(df_sub))
        
        ax.barh(y_pos - 0.2, df_sub["prop_a"], 0.4, label=sample_a, color="red", alpha=0.7)
        ax.barh(y_pos + 0.2, df_sub["prop_b"], 0.4, label=sample_b, color="blue", alpha=0.7)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"Subcluster {s}" for s in df_sub["subcluster"]], fontsize=9)
        ax.set_xlabel("Proportion within subcluster", fontsize=10)
        ax.set_ylabel("Subcluster", fontsize=10)
        ax.set_title(f"{ct}: Subcluster Composition\nSimilarity = {similarity:.3f}", fontsize=12, fontweight="bold")
        ax.legend(loc="upper right")
        ax.axvline(x=0.5, color="black", linestyle="--", alpha=0.5)
        ax.grid(axis="x", alpha=0.3)
        
        plt.tight_layout()
        out_file = os.path.join("figures", f"{args.prefix}_{ct}_subcluster_composition.png")
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        print(f"  Saved: {out_file}")
        plt.close()
    
    # Create similarity barplot only (no heatmap, no added numbers)
    valid_scores = {ct: score for ct, score in similarity_scores.items() if not np.isnan(score)}
    
    if valid_scores:
        celltype_list = list(valid_scores.keys())
        similarity_values = list(valid_scores.values())
        
        fig, ax = plt.subplots(figsize=(12, max(6, len(celltype_list) * 0.4)))
        
        colors = plt.cm.RdBu_r(np.array(similarity_values))
        ax.barh(celltype_list, similarity_values, color=colors)
        
        ax.set_xlim(0, 1)
        ax.set_xlabel('Similarity Score (1 = identical, 0 = completely different)', fontsize=10)
        ax.set_title(f'Subtype Composition Similarity Between {sample_a} and {sample_b}', fontsize=12, fontweight='bold')
        ax.axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        barplot_file = os.path.join("figures", f"{args.prefix}_subcluster_similarity_barplot.png")
        plt.savefig(barplot_file, dpi=300, bbox_inches='tight')
        print(f"\nSaved: {barplot_file}")
        plt.close()
    
    # Save results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(f"{args.prefix}_subcluster_composition.csv", index=False)
        print(f"\nSaved: {args.prefix}_subcluster_composition.csv")
        
        # Print summary
        print("\n" + "="*70)
        print("LOCAL STRUCTURE SIMILARITY SUMMARY")
        print("="*70)
        for ct in sorted(celltype_list, key=lambda x: similarity_scores[x], reverse=True):
            score = similarity_scores[ct]
            if score >= 0.8:
                rating = "Very similar"
            elif score >= 0.6:
                rating = "Moderately similar"
            elif score >= 0.4:
                rating = "Somewhat different"
            else:
                rating = "Very different"
            print(f"{ct:<25} {score:<10.3f} {rating}")
    
    print(f"\nDone")

if __name__ == "__main__":
    main()
