import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', required=True, help='DGE CSV file')
    parser.add_argument('--organism', default='mouse', help='Organism (default: mouse)')
    parser.add_argument('--resources', default='go', help='Gene set resource (default: go)')
    parser.add_argument('--prefix', default='GSEA_results', help='Prefix for output files')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    # Load data
    df = pd.read_csv(args.csv)
    print(f"Loaded {len(df)} genes")
    
    # Create ranking metric
    if 'pvals_adj' in df.columns:
        pvals = df['pvals_adj'].copy()
        min_p = pvals[pvals > 0].min() if any(pvals > 0) else 1e-300
        pvals = pvals.replace(0, min_p)
        neg_log_p = -np.log10(pvals)
        max_neg_log_p = neg_log_p.max()
        df['ranking'] = np.sign(df['logfoldchanges']) * (abs(df['logfoldchanges']) + neg_log_p / max_neg_log_p)
    else:
        df['ranking'] = df['logfoldchanges']
    
    # Prepare ranking
    ranking = df[['gene', 'ranking']].dropna().sort_values('ranking', ascending=False)
    ranking.columns = ['Gene', 'Score']
    
    dup_count = ranking['Score'].duplicated().sum()
    print(f"Duplicate scores: {dup_count} ({dup_count/len(ranking)*100:.1f}%)")
    
    if dup_count > 0:
        np.random.seed(42)
        ranking['Score'] = ranking['Score'] + np.random.normal(0, 1e-10, len(ranking))
    
    # Get GO gene sets for mouse
    from gseapy.parser import get_library_name
    libraries = get_library_name(organism=args.organism.capitalize())
    go_libs = [lib for lib in libraries if 'GO_Biological_Process' in lib]
    gene_set = go_libs[0] if go_libs else 'GO_Biological_Process_2023'
    print(f"Using gene set: {gene_set}")
    
    # Run GSEA
    results = gp.prerank(
        rnk=ranking,
        gene_sets=gene_set,
        permutation_num=1000,
        outdir=f'{args.prefix}_gsea_output',
        min_size=15,
        max_size=500,
        seed=42
    )
    
    # Convert numeric columns
    results_df = results.res2d
    numeric_cols = ['ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val', 'Tag %', 'Gene %']
    for col in numeric_cols:
        if col in results_df.columns:
            results_df[col] = pd.to_numeric(results_df[col], errors='coerce')
    
    results_df.to_csv(f'{args.prefix}_results.csv', index=False)
    print(f"\nResults saved to {args.prefix}_results.csv")
    print(f"Found {len(results_df)} pathways")
    
    # Clean NES column
    results_df = results_df.dropna(subset=['NES'])
    results_df = results_df.sort_values('NES', ascending=False)
    
    # Create dotplot
    print("\nGenerating dotplot...")
    
    # Get top 15 enriched and top 15 depleted
    top_enriched = results_df[results_df['NES'] > 0].head(15)
    top_depleted = results_df[results_df['NES'] < 0].tail(15)
    
    # Combine for plotting
    plot_df = pd.concat([top_enriched, top_depleted])
    
    # Add -log10(FDR) for color
    if 'FDR q-val' in plot_df.columns:
        plot_df['-log10(FDR)'] = -np.log10(plot_df['FDR q-val'] + 1e-10)
    else:
        plot_df['-log10(FDR)'] = 1
    
    # Create the dotplot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Normalize NES for dot size (size range 50-500)
    nes_range = plot_df['NES'].abs().max()
    sizes = 50 + (plot_df['NES'].abs() / nes_range) * 450
    
    # Color by -log10(FDR)
    scatter = ax.scatter(
        plot_df['NES'],
        range(len(plot_df)),
        s=sizes,
        c=plot_df['-log10(FDR)'],
        cmap='RdYlBu_r',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(FDR)', fontsize=10)
    
    # Add horizontal line at 0
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Set y-axis labels
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels([t[:50] + '...' if len(t) > 50 else t for t in plot_df['Term'].values], fontsize=8)
    
    # Labels and title
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_ylabel('Pathways', fontsize=12)
    ax.set_title(f'GSEA Results: Top 30 Pathways\n{args.organism.upper()} {args.resources.upper()}', 
                 fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='x')
    
    # Invert y-axis to show top enriched at top
    ax.invert_yaxis()
    
    # Add legend for dot sizes
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='|NES| = 1',
               markerfacecolor='gray', markersize=8, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='|NES| = 2',
               markerfacecolor='gray', markersize=13, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='|NES| = 3',
               markerfacecolor='gray', markersize=18, markeredgecolor='black'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', title='Dot size = |NES|', fontsize=8)
    
    plt.tight_layout()
    
    # Save figure
    output_path = 'figures/' + args.prefix + '_dotplot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.show()
    
    # Print top results to console
    print("\n" + "="*60)
    print("TOP 10 ENRICHED PATHWAYS (by NES):")
    print("="*60)
    for idx, row in results_df.head(10).iterrows():
        fdr = row['FDR q-val'] if 'FDR q-val' in results_df.columns else 'N/A'
        print(f"{row['Term'][:70]}: NES={row['NES']:.3f}, FDR={fdr}")
    
    print("\n" + "="*60)
    print("TOP 10 DEPLETED PATHWAYS (by NES):")
    print("="*60)
    for idx, row in results_df.tail(10).iterrows():
        fdr = row['FDR q-val'] if 'FDR q-val' in results_df.columns else 'N/A'
        print(f"{row['Term'][:70]}: NES={row['NES']:.3f}, FDR={fdr}")

if __name__ == "__main__":
    main()
