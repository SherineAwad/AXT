#!/usr/bin/env python3

import argparse
import scanpy as sc

marker_genes = {
    "Fibroblast": ["Prrx1", "Pdgfra", "Col1a1", "Col1a2", "Dcn", "Pi16"],
    "Endothelial": ["Cdh5", "Pecam1", "Kdr", "Emcn", "Erg"],
    "Macrophage": ["Adgre1", "Csf1r", "Cd68", "Mrc1"],
    "Keratinocyte": ["Krt14", "Krt5", "Krt17", "Krt6a"],
    "Osteoblast": ["Sp7", "Bglap", "Alpl", "Ibsp"],
    "Pericyte": ["Cspg4", "Pdgfrb", "Kcnj8", "Abcc9"],
    "SMC": ["Myh11", "Tagln", "Acta2", "Cnn1"],
    "Chondrocyte": ["Col2a1", "Acan", "Sox9"],
    "Schwann": ["Mbp", "Plp1", "Sox10", "S100b"],
    "T-cell": ["Cd3d", "Cd3e", "Cd8a", "Cd4"],
    "Osteoclast": ["Nfatc1", "Oscar", "Ctsk", "Acp5", "Mmp9", "Dcstamp"],
    "Synoviocyte": ["Prg4", "Ucma", "Gdf5", "Cilp2", "Frzb"],
    "Neutrophil": ["Ly6g", "S100a8", "S100a9", "Mpo"],
    "Lymphatic_Endothelial": ["Prox1", "Lyve1", "Pdpn", "Flt4"],
    "B-cell": ["Cd79a", "Cd19", "Ms4a1", "Ighm", "Pax5"],
    "Rspo3_Col23a1": ["Rspo3", "Col23a1"],
    "MSC": ["Lepr", "Cxcl12", "Pdgfra", "Eng"]
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# Collect all unique genes
all_genes = []
for genes in marker_genes.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))

# Filter and warn about missing genes
valid_genes = []
for gene in all_genes:
    if gene in adata.var_names:
        valid_genes.append(gene)
    else:
        print(f"WARNING: {gene} not found in dataset - skipping")

# Feature plot for each gene (gene name only in filename)
for gene in valid_genes:
    sc.pl.umap(adata, color=gene, save=f"_{args.prefix}_{gene}.png")

# Dotplot (uses full marker_genes dict with celltypes)
filtered_marker_genes = {}
for celltype, genes in marker_genes.items():
    filtered_genes = [g for g in genes if g in adata.var_names]
    if filtered_genes:
        filtered_marker_genes[celltype] = filtered_genes

sc.pl.dotplot(adata, filtered_marker_genes, groupby="leiden", save=f"_{args.prefix}_dotplot.png")
