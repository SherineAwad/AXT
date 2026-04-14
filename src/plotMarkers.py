#!/usr/bin/env python3

import argparse
import scanpy as sc

marker_genes = {
    "Fibroblast": ["Col1a1", "Col1a2", "Dcn", "Pdgfra", "Prrx1", "Pi16"],
    "Endothelial": ["Cdh5", "Pecam1", "Kdr", "Emcn", "Erg"],
    "Macrophage": ["Cd68", "Csf1r", "Adgre1", "Cd163", "Mrc1"],
    "Keratinocyte": ["Krt5", "Krt14", "Krt6a", "Krt17"],
    "Osteoblast": ["Sp7", "Runx2", "Alpl", "Bglap", "Ibsp"],
    "Pericyte": ["Pdgfrb", "Cspg4", "Rgs5", "Kcnj8", "Abcc9"],
    "SMC": ["Acta2", "Tagln", "Myh11", "Cnn1", "Acta1"],
    "Chondrocyte": ["Sox9", "Col2a1", "Acan", "Col10a1"],
    "Schwann": ["Sox10", "S100b", "Mbp", "Plp1", "Pmp22"],
    "T-cell": ["Cd3e", "Cd3d", "Cd4", "Cd8a", "Il2rb"],
    "Osteoclast": ["Ctsk", "Acp5", "Mmp9", "Dcstamp", "Traf6"],
    "Synoviocyte": ["Prg4", "Ucma", "Gdf5", "Cilp2", "Frzb"],
    "Neutrophil": ["Ly6g", "S100a8", "S100a9", "Cxcr2", "Mpo"],
    "Lymphatic_Endothelial": ["Prox1", "Lyve1", "Pdpn", "Flt4"],
    "B-cell": ["Cd19", "Ms4a1", "Cd79a", "Ighm", "Pax5"],
    "Rspo3_Col23a1": ["Rspo3", "Col23a1"],
    "MSC": ["Lepr", "Prrx1", "Cxcl12", "Pdgfra", "Eng", "Thy1"]
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
