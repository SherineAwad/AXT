# Project Overview: AXT Transplantation (2 wpa) scRNA-seq Analysis

This project studies tissue response **2 weeks post amputation (2 wpa)** using a transplantation model inspired by the axolotl (a highly regenerative animal).

Samples come from mouse and include two conditions:

- **Regenerative** → tissue attempting to regrow  
- **Non-regenerative** → normal healing (scar formation)

#### Objective
Compare:
- **Regeneration (regrowth)**  
- **Non-regeneration (failed healing)**

#### Key Goal
Identify **cell types and gene programs** that enable regeneration vs. lead to scarring.

## Filtering 

Cells are filtered to remove **low-quality cells**, **dead cells**, and **doublets**.

1. `n_genes_by_counts`
**Definition:** Number of unique genes detected per cell.

**Why filter:**
- Low → dead cells / empty droplets
- High → potential doublets

2. `total_counts`
**Definition:** Total RNA counts (UMIs) per cell.

**Why filter:**
- Low → poor RNA capture (noisy data)
- High → potential doublets

3. `pct_counts_mt`
**Definition:** Percentage of counts from mitochondrial genes.

**Why filter:**
- High → stressed or dying cells

### Pre-filtering 

![](figures/violin_axt_preQC.png?v=4)

### Post filtering 

![](figures/violin_axt_AfterQC.png?v=4) 


## Doublet Detection using Scrublet 


![](figures/axt_scrublet_scores.png?v=2) 

- Detected 82 doublets (0.3%) and were removed. 


## Analysing 

1. Normalization
`sc.pp.normalize_total(target_sum=1e4)`
- Makes all cells comparable by scaling total counts
- Each cell gets the same total (10,000 counts)

2. Log Transformation
`sc.pp.log1p()`
- Reduces effect of very large values
- Makes data more balanced and easier to analyze

3. Scaling
`sc.pp.scale(max_value=10)`
- Centers genes (mean = 0) and standardizes variance
- Clips extreme values to avoid outliers dominating

4. PCA (Dimensionality Reduction)
`sc.tl.pca()`
- Compresses data into main patterns (principal components)
- Keeps most important biological variation

5. Neighbors Graph
`sc.pp.neighbors()`
- Finds similar cells based on PCA
- Builds a graph of cell relationships

### 6. UMAP (Visualization)
`sc.tl.umap()`
- Projects cells into 2D space
- Similar cells cluster together visually

#### UMAP
![](figures/umap_axt_umap.png?v=5) 

#### Per sample UMAP 

<img src="figures/umap_axt_Reg.png?v=4" width="600" /> <img src="figures/umap_axt_nonReg.png?v=4" width="600" />


## Top PC genes 

#### PC1 positive genes
Tmsb4x, Srgn, Lcp1, Plek, Cd53, Ptprc, Laptm5, Cd52, Tyrobp, Mir142hg

#### PC1 negative genes
Nedd4, Rpl39l, EGFP, AI506816, Tpm1, Cald1, Serpinh1, Npm1, Cks1b, Col6a3

#### PC2 positive genes
Igfbp7, Zbtb20, Serping1, Plpp3, Plac9, Ebf1, Dcn, Ccdc80, Fstl1, Lama2

#### PC2 negative genes
H2az1, Pfn1, Fau, Rpl27, Rpl38, Rpl35a, Fabp5, Actb, Rps29, Rpl34

#### PC3 positive genes
Tpt1, Itm2b, S100a6, Eif1, Gabarap, Ubb, Gpx1, Myl6, Ppib, Ftl1

#### PC3 negative genes
Taco1, Cdk8, Kcnq5, Gphn, Cmss1, Snd1, Lars2, Camk1d, Gm12610, Smyd3
`
## Clustering 

`sc.tl.leiden(resolution=1.0)`

#### What it does
- Groups cells into **clusters** based on similarity
- Uses the **neighbors graph** (cells that are close to each other)

#### How to think about it
- Cells with similar gene expression → same cluster  
- Each cluster ≈ **cell type or cell state**

#### Resolution parameter
- **Low resolution (e.g. 0.3)** → fewer, larger clusters  
- **High resolution (e.g. 1.0+)** → more, smaller clusters  

## Clustering (Leiden)

`sc.tl.leiden(resolution=1.0)`

#### What it does
- Groups cells into **clusters** based on similarity
- Uses the **neighbors graph** (cells that are close to each other)

#### How to think about it
- Cells with similar gene expression → same cluster  
- Each cluster ≈ **cell type or cell state**

#### Resolution parameter
- **Low resolution (e.g. 0.3)** → fewer, larger clusters  
- **High resolution (e.g. 1.0+)** → more, smaller clusters  


![](figures/umap_axt_leiden.png?v=4) 


## QC per Leiden Cluster

QC metrics were visualized across Leiden clusters to assess cluster quality.
Check QC per cluster to spot and remove **low-quality or suspicious cell groups**

<img src="figures/violin_axt_QC_n_genes_by_counts.png?v=5" width="33%" /><img src="figures/violin_axt_QC_total_counts.png?v=5" width="33%" /> <img src="figures/violin_axt_QC_pct_counts_mt.png?v=5" width="33%" />

## Marker genes 

To annotate cell types, I **manually selected (guessed) canonical marker genes** based on known biology from literature.

These markers were then used to validate and assign identities to clusters.


## Marker Gene Sets

```python
marker_genes = {
    "Fibroblast": ["Prrx1", "Pdgfra", "Col1a1", "Dcn", "Pi16", "Cd34"],
    "Endothelial": ["Esam", "Flt1", "Vwf", "Plvap","Cdh5", "Pecam1", "Kdr", "Emcn", "Erg", "Cd34"],
    "Macrophage": ["Adgre1", "Csf1r", "Cd68", "Mrc1", "Cd163"],
    "Keratinocyte": ["Krt14", "Krt5", "Epcam", "Cdh1", "Krt17", "Dsg3"],
    "Osteoblast": ["Runx2", "Postn", "Mmp13", "Spp1", "Dmp1","Sp7", "Bglap", "Alpl", "Ibsp", "Col1a1"],
    "Pericyte": ["Cspg4", "Pdgfrb", "Kcnj8", "Abcc9", "Rgs5"],
    "SMC": ["Myh11", "Tagln", "Acta2", "Cnn1", "Des"],
    "Chondrocyte": ["Col2a1", "Acan", "Sox9", "Matn1", "Col9a1", "Comp"],
    "Schwann": ["Mbp", "Plp1", "Sox10", "S100b", "Pmp22"],
    "T-cell": ["Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd28"],
    "Osteoclast": ["Ctsk", "Acp5", "Calcr", "Oscar", "Nfatc1", "Dcstamp", "Tnfrsf11a"],
    "Synoviocyte": ["Prg4", "Ucma", "Gdf5", "Cilp2", "Frzb"],
    "Neutrophil": ["Ly6g", "S100a8", "S100a9", "Mpo", "Csf3r"],
    "Lymphatic_Endothelial": ["Prox1", "Lyve1", "Pdpn", "Flt4", "Pecam1"],
    "B-cell": ["Cd79a", "Cd19", "Ms4a1", "Ighm", "Pax5", "Cd22"],
    "Rspo3_Col23a1": ["Rspo3", "Col23a1"],
    "MSC": ["Lepr", "Cxcl12", "Ngfr", "Nes", "Cd44", "Scf"],
    "Osteosarcoma": ["EGFP"],
    "Nail_Epithelium": ["Lgr6", "Sp6", "Sp8"],
    "Sweat_glands": ["Aqp5","Scnn1a","Scnn1b","Scnn1g","Krt19","Krt7","Krt8","Krt18","Krt5","Krt14","Foxa1"],
    "Mast cell": ["Kit", "Cpa3", "Tpsab1", "Tpsb2", "Ms4a2", "Hdc", "Hpgds", "Mcpt8", "Cd200r3", "Ccr3"]
}

```

### Dotplot for marker genes 

![](figures/dotplot__axt_dotplot.png?v=7) 


### Feature plots for marker genes 

<img src="figures/umap_axt_Rgs5.png?v=4" width="33%" /><img src="figures/umap_axt_Krt5.png?v=4" width="33%" /><img src="figures/umap_axt_Ms4a1.png?v=4" width="33%" />
<img src="figures/umap_axt_Prrx1.png?v=4" width="33%" /><img src="figures/umap_axt_Krt14.png?v=4" width="33%" /><img src="figures/umap_axt_Krt17.png?v=4" width="33%" />
<img src="figures/umap_axt_Sp6.png?v=4" width="33%" /><img src="figures/umap_axt_Nfatc1.png?v=4" width="33%" /><img src="figures/umap_axt_Cd68.png?v=4" width="33%" />
<img src="figures/umap_axt_Lgr6.png?v=4" width="33%" /><img src="figures/umap_axt_Cd22.png?v=4" width="33%" /><img src="figures/umap_axt_Ibsp.png?v=4" width="33%" />
<img src="figures/umap_axt_Mbp.png?v=4" width="33%" /><img src="figures/umap_axt_Cd79a.png?v=4" width="33%" /><img src="figures/umap_axt_Cd8a.png?v=4" width="33%" />
<img src="figures/umap_axt_Dcn.png?v=4" width="33%" /><img src="figures/umap_axt_Rspo3.png?v=4" width="33%" /><img src="figures/umap_axt_Cdh1.png?v=4" width="33%" />
<img src="figures/umap_axt_Oscar.png?v=4" width="33%" /><img src="figures/umap_axt_Cd44.png?v=4" width="33%" /><img src="figures/umap_axt_Mrc1.png?v=4" width="33%" />
<img src="figures/umap_axt_Col2a1.png?v=4" width="33%" /><img src="figures/umap_axt_Ucma.png?v=4" width="33%" /><img src="figures/umap_axt_Cd34.png?v=4" width="33%" />
<img src="figures/umap_axt_Ngfr.png?v=4" width="33%" /><img src="figures/umap_axt_Comp.png?v=4" width="33%" /><img src="figures/umap_axt_Epcam.png?v=4" width="33%" />
<img src="figures/umap_axt_Pax5.png?v=4" width="33%" /><img src="figures/umap_axt_Ctsk.png?v=4" width="33%" /><img src="figures/umap_axt_Lepr.png?v=4" width="33%" />
<img src="figures/umap_axt_Cd28.png?v=4" width="33%" /><img src="figures/umap_axt_Sox10.png?v=4" width="33%" /><img src="figures/umap_axt_Cspg4.png?v=4" width="33%" />
<img src="figures/umap_axt_Nes.png?v=4" width="33%" /><img src="figures/umap_axt_Cilp2.png?v=4" width="33%" /><img src="figures/umap_axt_Mpo.png?v=4" width="33%" />
<img src="figures/umap_axt_Prg4.png?v=4" width="33%" /><img src="figures/umap_axt_Tnfrsf11a.png?v=4" width="33%" /><img src="figures/umap_axt_Lyve1.png?v=4" width="33%" />
<img src="figures/umap_axt_Krt19.png?v=4" width="33%" /><img src="figures/umap_axt_Pecam1.png?v=4" width="33%" /><img src="figures/umap_axt_Pi16.png?v=4" width="33%" />
<img src="figures/umap_axt_Kcnj8.png?v=4" width="33%" /><img src="figures/umap_axt_Scnn1b.png?v=4" width="33%" /><img src="figures/umap_axt_Pdgfrb.png?v=4" width="33%" />
<img src="figures/umap_axt_Acp5.png?v=4" width="33%" /><img src="figures/umap_axt_Dcstamp.png?v=4" width="33%" /><img src="figures/umap_axt_Sp8.png?v=4" width="33%" />
<img src="figures/umap_axt_Acta2.png?v=4" width="33%" /><img src="figures/umap_axt_Frzb.png?v=4" width="33%" /><img src="figures/umap_axt_Cd163.png?v=4" width="33%" />
<img src="figures/umap_axt_Cd19.png?v=4" width="33%" /><img src="figures/umap_axt_Pdgfra.png?v=4" width="33%" /><img src="figures/umap_axt_S100a8.png?v=4" width="33%" />
<img src="figures/umap_axt_Cnn1.png?v=4" width="33%" /><img src="figures/umap_axt_Krt7.png?v=4" width="33%" /><img src="figures/umap_axt_Plp1.png?v=4" width="33%" />
<img src="figures/umap_axt_Myh11.png?v=4" width="33%" /><img src="figures/umap_axt_Prox1.png?v=4" width="33%" /><img src="figures/umap_axt_Gdf5.png?v=4" width="33%" />
<img src="figures/umap_axt_Sox9.png?v=4" width="33%" /><img src="figures/umap_axt_Acan.png?v=4" width="33%" /><img src="figures/umap_axt_Col9a1.png?v=4" width="33%" />
<img src="figures/umap_axt_Scnn1g.png?v=4" width="33%" /><img src="figures/umap_axt_Pmp22.png?v=4" width="33%" /><img src="figures/umap_axt_Cdh5.png?v=4" width="33%" />
<img src="figures/umap_axt_Scnn1a.png?v=4" width="33%" /><img src="figures/umap_axt_Ly6g.png?v=4" width="33%" /><img src="figures/umap_axt_S100a9.png?v=4" width="33%" />
<img src="figures/umap_axt_Flt4.png?v=4" width="33%" /><img src="figures/umap_axt_Foxa1.png?v=4" width="33%" /><img src="figures/umap_axt_Pdpn.png?v=4" width="33%" />
<img src="figures/umap_axt_Adgre1.png?v=4" width="33%" /><img src="figures/umap_axt_Des.png?v=4" width="33%" /><img src="figures/umap_axt_Calcr.png?v=4" width="33%" />
<img src="figures/umap_axt_Cxcl12.png?v=4" width="33%" /><img src="figures/umap_axt_Abcc9.png?v=4" width="33%" /><img src="figures/umap_axt_Cd3d.png?v=4" width="33%" />
<img src="figures/umap_axt_Matn1.png?v=4" width="33%" /><img src="figures/umap_axt_Col23a1.png?v=4" width="33%" /><img src="figures/umap_axt_Ighm.png?v=4" width="33%" />
<img src="figures/umap_axt_Dsg3.png?v=4" width="33%" /><img src="figures/umap_axt_Alpl.png?v=4" width="33%" /><img src="figures/umap_axt_Bglap.png?v=4" width="33%" />
<img src="figures/umap_axt_Cd3e.png?v=4" width="33%" /><img src="figures/umap_axt_Erg.png?v=4" width="33%" /><img src="figures/umap_axt_S100b.png?v=4" width="33%" />
<img src="figures/umap_axt_Csf3r.png?v=4" width="33%" /><img src="figures/umap_axt_Aqp5.png?v=4" width="33%" /><img src="figures/umap_axt_Csf1r.png?v=4" width="33%" />
<img src="figures/umap_axt_Tagln.png?v=4" width="33%" /><img src="figures/umap_axt_Emcn.png?v=4" width="33%" /><img src="figures/umap_axt_Krt18.png?v=4" width="33%" />
<img src="figures/umap_axt_Sp7.png?v=4" width="33%" /><img src="figures/umap_axt_Kdr.png?v=4" width="33%" /><img src="figures/umap_axt_Krt8.png?v=4" width="33%" />
<img src="figures/umap_axt_Rgs5.png?v=4" width="33%" /><img src="figures/umap_axt_Krt5.png?v=4" width="33%" /><img src="figures/umap_axt_Ms4a1.png?v=4" width="33%" />
<img src="figures/umap_axt_Prrx1.png?v=4" width="33%" /><img src="figures/umap_axt_Krt14.png?v=4" width="33%" /><img src="figures/umap_axt_Krt17.png?v=4" width="33%" />
<img src="figures/umap_axt_Kit.png?v=4" width="33%" /><img src="figures/umap_axt_Cpa3.png?v=4" width="33%" /><img src="figures/umap_axt_Tpsab1.png?v=4" width="33%" />
<img src="figures/umap_axt_Tpsb2.png?v=4" width="33%" /><img src="figures/umap_axt_Ms4a2.png?v=4" width="33%" /><img src="figures/umap_axt_Hdc.png?v=4" width="33%" />
<img src="figures/umap_axt_Hpgds.png?v=4" width="33%" /><img src="figures/umap_axt_Mcpt8.png?v=4" width="33%" /><img src="figures/umap_axt_Cd200r3.png?v=4" width="33%" />
<img src="figures/umap_axt_Ccr3.png?v=4" width="33%" />

`
## Annotations 

We used the above marker genes to annotate our celltypes as follows:

![](figures/umap_axt_celltype.png?v=8)

![](figures/umap_axt_celltypeON.png?v=8)

## Some stats for samples and cells 

| Celltype                     | Reg   | nonReg | Total |
|------------------------------|------:|-------:|------:|
| B-Cells                      |   259 |    291 |   550 |
| Basophil                     |    30 |    109 |   139 |
| Chondrocyte                  |   185 |      3 |   188 |
| Endothelial                  |  1513 |   1076 |  2589 |
| Fibroblast                   |  2759 |    827 |  3586 |
| Lymphatic_Endothelial        |   147 |    203 |   350 |
| Macrophage                   |  1774 |   3221 |  4995 |
| NailEpithelium_Keratinocyte  |   663 |    491 |  1154 |
| Neutrophil                   |   221 |    424 |   645 |
| OsteoProgenitor              |  1341 |   1797 |  3138 |
| Osteoblast                   |    80 |     39 |   119 |
| Osteoclast                   |   136 |    561 |   697 |
| Osteosarcoma                 |  1139 |   9034 | 10173 |
| SMC_Pericyte                 |   553 |    540 |  1093 |
| Schwann                      |   176 |     54 |   230 |
| SweatGland                   |   215 |      2 |   217 |
| T-Cells                      |   523 |   1014 |  1537 |
| **Total**                    | 11714 |  19686 | 31400 |


![](figures/axt_cell_ratios.png?v=9) 

## Decontaminate Osteosarcoma 

#### Cdkn2a > 2.0 and EGFP < 1.0 
We decontaminate the Osteosarcoma population by removing cells that express Cdkn2a at levels > 2.0, provided that their EGFP expression is < 1.0.

![](figures/umap_axt_before_minExp2.0_maxEGFP1.0.png?v=1)

![](figures/umap_axt_after_minExp2.0_maxEGFP1.0.png?v=1)

| Stage   | Cells  |
|---------|-------:|
| Before  | 10,173 |
| After   | 9,158  |
| Removed | 1,015  |


#### Cdkn2a > 2.0 and EGFP < 0.5 
We decontaminate the Osteosarcoma population by removing cells that express Cdkn2a at levels > 2.0, provided that their EGFP expression is < 0.5.

![](figures/umap_axt_before_minExp2.0_maxEGFP0.5.png?v=1) 

![](figures/umap_axt_after_minExp2.0_maxEGFP0.5.png?v=1)

| Metric   | Cells |
|----------|------:|
| Before   | 10,173 |
| After    | 9,600  |
| Removed  | 573    |



#### Cdkn2a > 1.5 and EGPF < 0.5 

We decontaminate the Osteosarcoma population by removing cells that express Cdkn2a at levels > 1.5, provided that their EGFP expression is < 0.5.

![](figures/umap_axt_before_minExp1.5_maxEGFP0.5.png?v=1)

![](figures/umap_axt_after_minExp1.5_maxEGFP0.5.png?v=1)

| Metric   | Cells |
|----------|------:|
| Before   | 10,173 |
| After    | 9,417  |
| Removed  | 756    |


#### Cdkn2a > 1.0 and EGFP < 0.5 

We decontaminate the Osteosarcoma population by removing cells that express Cdkn2a at levels > 1.0, provided that their EGFP expression is < 0.5.

![](figures/umap_axt_before_minExp1.0_maxEGFP0.5.png?v=1)

![](figures/umap_axt_after_minExp1.0_maxEGFP0.5.png?v=1)

| Metric   | Cells |
|----------|------:|
| Before   | 10,173 |
| After    | 9,259  |
| Removed  | 914    |


## Proceeding decontamination with CDkn2a > 1.0 and EGFP < 0.5 


### Cellration update after cells removal

| Celltype                 | Reg   | nonReg | Total |
|--------------------------|------:|------:|------:|
| B-Cells                  | 259   | 291   | 550   |
| Basophil                 | 30    | 109   | 139   |
| Chondrocyte              | 185   | 3     | 188   |
| Endothelial              | 1513  | 1076  | 2589  |
| Fibroblast               | 2759  | 827   | 3586  |
| Lymphatic_Endothelial    | 147   | 203   | 350   |
| Macrophage               | 1774  | 3221  | 4995  |
| NailEpithelium_Keratinocyte | 663 | 491   | 1154  |
| Neutrophil               | 221   | 424   | 645   |
| OsteoProgenitor          | 1341  | 1797  | 3138  |
| Osteoblast               | 80    | 39    | 119   |
| Osteoclast               | 136   | 561   | 697   |
| Osteosarcoma             | 1033  | 8226  | 9259  |
| SMC_Pericyte             | 553   | 540   | 1093  |
| Schwann                  | 176   | 54    | 230   |
| SweatGland               | 215   | 2     | 217   |
| T-Cells                  | 523   | 1014  | 1537  |
| **Total**                | 11608 | 18878 | 30486 |

![](figures/axt_decontaminated_cell_ratios.png?v=1)

## 🚨⚠️ **We have to repeat the downstream analysis after the deconatimantion step** ⚠️  🚨
`

## Zooming on Osteosarcoma 

- We subset Osteosarcoma, below is a summary: 

| Metric  | Value |
|---------|------:|
| Cells   | 9,259 |
| Genes   | 27,808 |

| Sample  | Count |
|---------|------:|
| nonReg  | 8,226 |
| Reg     | 1,033 |

![](figures/Osteosarcoma_umap.png?v=2)

![](figures/Osteosarcoma_perSample_umap.png?v=2)


### Differential gene expression in Osteosarcoma 

Differential gene expression (DGE) was performed using a Wilcoxon rank-sum test on log-normalised expression data to compare sample groups against a defined reference condition (nonReg). Genes were filtered based on an adjusted p-value threshold of 0.05 to retain statistically significant changes. For heatmap visualisation, the top 100 genes per direction (up- and downregulated) were selected based on adjusted p-value ranking with log fold-change used as a secondary criterion. For the volcano plot, all genes were displayed, while the top 10 most statistically significant genes per direction were highlighted and labelled. Gene ranking for selection in both visualisations was primarily driven by adjusted p-value, with log fold-change used to refine ordering within directionally separated groups.

#### Top 50 genes 
![](figures/Osteosarcoma_heatmap.png?v=6) 

#### Volcano plot higlighting genes of interest 

![](figures/Osteosarcoma_volcano.png?v=8) 

#### Dotplot for genes of interest in Osteosarcoma 

![](figures/dotplot__Osteosarcoma_dotplot.png?v=1)

#### Violin plots for genes of interest in Osteosarcoma 

<img src="figures/Osteosarcoma_Sparc_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Omd_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Cdkn1a_violin.png?v=2" width="33%" />
<img src="figures/Osteosarcoma_Alpl_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Wif1_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Panx3_violin.png?v=2" width="33%" />
<img src="figures/Osteosarcoma_Kazald1_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Dcn_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Nbl1_violin.png?v=2" width="33%" />
<img src="figures/Osteosarcoma_Bglap_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Bglap2_violin.png?v=2" width="33%" /><img src="figures/Osteosarcoma_Aspn_violin.png?v=2" width="33%" />
<img src="figures/Osteosarcoma_Notum_violin.png?v=2" width="33%" />

#### Feature plots for genes of interest in Osteosarcoma

<img src="figures/Osteosarcoma_Notum_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Bglap2_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Aspn_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Bglap_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Nbl1_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Panx3_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Kazald1_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Dcn_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Wif1_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Cdkn1a_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Alpl_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Omd_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Sparc_featureplot.png?v=1" width="33%" />

[📊📊 Download here DGE for Osteosarcoma Reg vs non Reg with adj-pvalue \<0.05](https://docs.google.com/spreadsheets/d/1NbpSDUpz78nXY_4nH0Ksp0ip3VMUfFLsjfP-Aujh6GA/edit?usp=sharing)

### A look into Proliferation Genes

To rule out proliferation as a driver of the Reg vs nonReg transcriptional difference in Osteosarcoma, we examined expression of canonical proliferation markers.

![](figures/dotplot__Proliferation_Osteosarcoma_dotplot.png?v=2)

### A look into Inflammation Genes 

To rule out inflammation as a driver of the Reg vs nonReg transcriptional difference in Osteosarcoma, we examined expression of canonical inflammatory markers.

![](figures/dotplot__Inflammation_Osteosarcoma_dotplot.png?v=1) 

## Pathways and GO Enrichment for  Osteosarcoma 


### GSEA GO Enrichments 

We calculate a ranking score for every gene using: sign(logFC) × (abs(logFC) + normalized p-value). This puts upregulated genes at the top and downregulated at the bottom, with logFC as the primary factor and p-values as tie-breakers. Genes with identical scores get tiny random noise to break ties.

For each GO pathway, we walk through the ranked gene list and asks: are genes from this pathway clustered near the top (upregulated) or bottom (downregulated)? It calculates an Enrichment Score (ES), then normalizes it to account for pathway size, producing the Normalized Enrichment Score (NES). Positive NES = pathway activated. Negative NES = pathway suppressed. The script does 1000 random permutations to calculate FDR q-values.

We plot the top 15 enriched (highest NES) and top 15 depleted (lowest NES) pathways. X-axis = NES (positive right, negative left). Dot size = |NES| (bigger = stronger). Dot color = -log10(FDR) (darker = more significant). The plot saves to figures/ directory.

![](figures/OsteosarcomaGSEA_dotplot.png?v=1)

[Download here the full list of Osteosarcoma GSEA](https://docs.google.com/spreadsheets/d/1bJ09-pqIcwq0liNOfoOlPIxs1TCgSPrBq1jRFR21isk/edit?usp=sharing)

### g:Profiler enrichment analysis 

g:Profiler is a tool for interpreting gene lists by mapping them onto known biological knowledge bases such as KEGG, Reactome, Gene Ontology, CORUM, and WikiPathways.

It answers a simple question:

> Are certain biological pathways or functions over-represented in a given list of genes?

Instead of focusing on individual genes, it summarizes them into higher-level biological processes.

Each pathway is evaluated based on how strongly it overlaps with the input gene list compared to what would be expected by chance. Pathways are then **ranked by enrichment strength**, so the most biologically relevant processes appear at the top.

This is **over-representation analysis (ORA)**:
it tests whether predefined biological categories are enriched in a selected gene set, providing a compact functional interpretation of differential expression results.

We performed over‑representation enrichment analysis using g:Profiler against the following functional annotation databases:

- **KEGG**: Manually curated pathway maps for metabolism, signaling, and diseases.
- **REAC (Reactome)**: Peer‑reviewed, hierarchical pathway database covering molecular processes.
- **GO (Gene Ontology)**: Three structured ontologies (BP, CC, MF) describing gene function.
- **WP (WikiPathways)**: Community‑curated, open platform for biological pathways.
- **CORUM**: Experimentally validated mammalian protein complexes


#### Using KEGG, REACTOME, and WP
![](figures/Osteosarcoma_KEGG,REAC,WP_enrichment.png?v=4)

[Download here the Osteosarcoma KEGG/REAC/WP enrichments](https://docs.google.com/spreadsheets/d/1Afs9llmco7xB1_ZJHWC2891UMkQCqmIniJyxU13Kp5s/edit?usp=sharing)

#### Using GO
![](figures/Osteosarcoma_GO_enrichment.png?v=4) 


#### A selective GO plot for interesting enrichments 

![](figures/axt_GO_enrichment_dotplot.png?v=1)

[Download here the Osteosarcoma GO enrichments](https://docs.google.com/spreadsheets/d/1WMmBSIl0G6ohwbj5UM__B_VqlxE4wHztsTT1ixJeZrQ/edit?usp=sharing)

#### Using CORUM 
![](figures/Osteosarcoma_CORUM_enrichment.png?v=4)  

[Download here the Osteosarcoma CORUM enrichments](https://docs.google.com/spreadsheets/d/1rYkPxymde6Q7bsLjTkHMDk3TlbbCveYegTtDbciTq_s/edit?usp=sharing)

#### How to read Go profiler results 

| Field | Meaning |
|-------|---------|
| **source** | The database or ontology source where this term comes from (like GO:BP for Gene Ontology Biological Process, KEGG for Kyoto Encyclopedia of Genes and Genomes, REAC for Reactome pathways) |
| **native** | The original identifier of the term from its source database (like a GO ID or KEGG pathway ID), which is the internal code used by the source database |
| **name** | The human-readable description of the term (for example, "apoptotic process" or "cell cycle"), which is the actual biological name you would use when writing about your results |
| **p_value** | The statistical significance score from the enrichment test, where lower values (closer to zero) mean the enrichment is more unlikely to have occurred by chance |
| **significant** | A boolean flag (True/False) indicating whether this term passed the significance threshold, based on the g:Profiler's internal multiple testing correction |
| **description** | A longer text explanation of what the term means biologically, which is similar to the "name" field but may contain additional details |
| **term_size** | The total number of genes annotated to this term in the entire genome/organism, representing the full size of this gene set in the background |
| **query_size** | The number of genes you submitted for enrichment analysis (your input gene list), which is the total genes you asked g:Profiler to test |
| **intersection_size** | The number of genes from your query that overlap with this term, meaning how many of your input genes are annotated to this specific pathway or process |
| **effective_domain_size** | The total number of genes considered as the statistical background after filtering, typically the number of genes that have at least one annotation in the sources you selected |
| **precision** | The proportion of your query genes in this term relative to the term size, essentially measuring how specific this term is to your gene set (precision = intersection_size / term_size) |
| **recall** | The proportion of your query genes in this term relative to your total query size, essentially measuring how well this term covers your input gene list (recall = intersection_size / query_size) |
| **query** | The specific subset of genes from your input that map to this term, usually shown as a comma-separated list of gene identifiers |
| **parents** | The parent terms in the ontology hierarchy, showing broader categories that contain this term (for example, a parent might be "biological regulation" for a child term like "cell cycle regulation") |
| **celltype** | Your custom added field identifying which cell type or experimental condition this enrichment came from, useful when analyzing multiple groups separately |
| **direction** | Your custom added field indicating whether this enrichment came from up-regulated genes (positive log fold change) or down-regulated genes (negative log fold change) |
| **gene_count** | Your custom added field that is exactly the same as intersection_size, just renamed for clarity in your plotting code |
| **neg_log10_pval** | Your custom calculated field which is -log10(p_value), where small p-values become larger positive numbers (p=0.001 becomes 3, p=0.000001 becomes 6), making visualization easier on plots |

## Zooming on Fibroblast 

We subset Fibroblast cell type, below is a summary:

| Metric  | Value |
|---------|------:|
| Cells   | 3,586 |
| Genes   | 27,808 |

| Sample  | Count |
|---------|------:|
| Reg     | 2,759 |
| nonReg  | 827 |

![](figures/Fibroblast_umap.png?v=2)

![](figures/Fibroblast_perSample_umap.png?v=2)

### Differential gene expression in Fibroblast

Similar to what we did with Osteosarcoma, Top N genes are selected by first sorting all genes by adjusted p-value in ascending order (most statistically significant first), with log fold-change used only as a secondary criterion to break ties while preserving directionality. The top-ranked genes are then taken globally from this ordered list.

#### Top 50 
![](figures/Fibroblast_heatmap.png?v=6)

Now overlapping with Secreted genes 
  
#### Top 10 genes only are labeled 
![](figures/Fibroblast_volcano.png?v=7)

#### Dotplot for genes of interest in Fibroblast  

![](figures/Fibroblast_dotplot.png?v=7)


#### Feature plots for genes of interest in Fibroblast

<img src="figures/Fibroblast_feature_Bmp2.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Aldh3a1.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Bmp5.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Bmp4.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Cdk8.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Cd81.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Camk1d.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Dpt.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Cmss1.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Foxp1.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Eif1.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Fxyd1.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Ftl1.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Il31ra.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Hexb.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Gphn.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Mfap4.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Lars2.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Tmem158.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Qpct.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_Vim.png?v=4" width="33%" />
<img src="figures/Fibroblast_feature_Ubb.png?v=4" width="33%" /><img src="figures/Fibroblast_feature_mt-Nd2.png?v=4" width="33%" />

#### Violin plots for genes of interest in Fibroblast 

<img src="figures/Fibroblast_violin_Bmp4.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Bmp2.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Aldh3a1.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Cdk8.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Cd81.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Camk1d.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Bmp5.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Eif1.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Dpt.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Cmss1.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Gphn.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Fxyd1.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Ftl1.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Foxp1.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Lars2.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Il31ra.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Hexb.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Tmem158.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Qpct.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Mfap4.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_Vim.png?v=4" width="33%" />
<img src="figures/Fibroblast_violin_Ubb.png?v=4" width="33%" /><img src="figures/Fibroblast_violin_mt-Nd2.png?v=4" width="33%" />


[📊📊Download Fibroblast DGE Reg vs non Reg with adj-pvalue \<0.05](https://docs.google.com/spreadsheets/d/11xskz0pwZ4XTII-RSRkLTfdnz2vGlNxauW287BD5dFs/edit?usp=sharing) 


### A look into Proliferation Genes

As we did in Osteosarcoma, to rule out proliferation as a driver of the Reg vs nonReg transcriptional difference in Fibroblast, we examined expression of canonical proliferation markers.

![](figures/dotplot__Proliferation_Fibroblast_dotplot.png?v=2)

### A look into Inflammation Genes 

As we did in Osteosarcoma, to rule out inflammation as a driver of the Reg vs nonReg transcriptional difference in Fibroblast, we examined expression of canonical inflammatory markers.

![](figures/dotplot__Inflammation_Fibroblast_dotplot.png?v=1)


### GOProfiler: Pathways and GO enrichments for Fibroblast  

We followed the same approach as in Osteosarcoma 

#### Using KEGG, REACTOME, and WP
![](figures/Fibroblast_KEGG,REAC,WP_enrichment.png?v=3)

[Download here the Fibroblast KEGG/REAC/WP enrichments](https://docs.google.com/spreadsheets/d/1kMVLqQM8qotiPLRciVcI_GueAI-u3T8Lw3sUHRsq15g/edit?usp=sharing)

#### Using GO
![](figures/Fibroblast_GO_enrichment.png?v=3) 

[Download here the Fibroblast GO enrichments](https://docs.google.com/spreadsheets/d/1Adn8g3giE3Ff7IKW9ME_RGSMpSarFGi7XMSW_ZYTjR4/edit?usp=sharing)

#### Using CORUM
![](figures/Fibroblast_CORUM_enrichment.png?v=3)   

[Download here the Fibroblast CORUM Enrichments](https://docs.google.com/spreadsheets/d/1TJwrrmCdz7Pft5gETs2HbkWPw-JPNQXV159pKsq5jcE/edit?usp=sharing)

[Go to how to read go profiler results section](#how-to-read-go-profiler-results)

## Now zooming on OsteoProgenitor 

We subset the OsteoProgenitor and a summary is below: 

| Metric | Value |
|--------|------:|
| Cells  | 3,138 |
| Genes  | 27,808 |

| Sample | Count |
|--------|------:|
| nonReg | 1,797 |
| Reg    | 1,341 |

![](figures/OsteoProgenitor_umap.png?v=2)

![](figures/OsteoProgenitor_perSample_umap.png?v=2)


### Differential gene expression in OsteoProgenitor 

We followed the same approach we did with Osteosarcoma (same script)

#### Top 50 genes
![](figures/OsteoProgenitor_heatmap.png?v=5)


## Global Regenerative vs non Regenerative
 
Here, we followed the same approach we did in Osteosarcoma for dge, but we compare Reg vs non Reg. 

#### Top 50 
![](figures/axt_heatmap.png?v=5)


[📊📊 Downdload here: Reg vs nonReg DGE with cutoff \<0.05](https://docs.google.com/spreadsheets/d/11-FSqCPBHs091N812u3h6CLEcuL40yKE7gKELVwnX84/edit?usp=sharing)


## Zooming into celltypes

We performed differential gene expression (DGE) analysis to identify genes differentially expressed across cell types. For each cell type, we applied a Wilcoxon rank-sum test using Scanpy’s rank_genes_groups function with celltype as the grouping variable, comparing each cell type against all remaining cells (one-vs-rest approach). We extracted log fold-changes and adjusted p-values for each comparison, and filtered results using an adjusted p-value threshold of 0.05 to retain statistically significant genes. For downstream visualisation, we selected the top 10 upregulated and top 10 downregulated genes per cell type (--N 10), ranking genes primarily by adjusted p-value with log fold-change used as a secondary criterion within each direction. The final results were visualised in a heatmap displaying log fold-changes of the selected genes across all cell types.

### Here, we used N=5 per celltype 

![](figures/axt_celltype_heatmap.png?v=11)

[📊📊 Downdload here: per celltype DGE with cutoff \<0.05](https://docs.google.com/spreadsheets/d/1DnR2domO66krErCqKF8BEU3zf5IzG5T6EZ1Ioekt9lE/edit?usp=sharing) 


## 🧬 Cell–Cell Interaction Analysis (LIANA)

This pipeline identifies potential **ligand–receptor interactions** from single-cell RNA-seq data using the LIANA framework. LIANA integrates multiple established methods into a single consensus scoring system.

#### Top lrScores interactions with minimum |logFC| 0.0 (Mostly rank by lrScore)

![](figures/axt_liana_dotplot.png?=7)

[📊📊Download full list of liana ligand-receptors interactions](https://docs.google.com/spreadsheets/d/1R2F82YcbmA1HZVSJf199n3f2_aNu_c1T5DaK8LR1iwk/edit?usp=sharing)

### Ligand–receptor inference (LIANA)
The script runs `liana.mt.rank_aggregate`, which integrates multiple tools:

- 📡 **CellPhoneDB** → permutation-based statistical testing of ligand–receptor pairs  
- 🧠 **CellChat** → probabilistic modeling of signaling pathway activity  
- 🔗 **NATMI** → network-based interaction scoring  
- 🌐 **Connectome** → expression-based interaction networks  
- 🧬 **SingleCellSignalR** → statistical inference of signaling relationships  

#### Consensus scoring
LIANA combines all methods into unified metrics:

- **`lrscore`** → Overall interaction confidence (0–1)  
  → Higher = stronger agreement across methods

- **`lr_logfc`** → Differential interaction strength  
  → Positive = enriched in condition A (e.g. Reg)  
  → Negative = enriched in condition B (e.g. nonReg)

#### Filtering and ranking
- Filter interactions by **|lr_logfc| threshold (absolute log fold-change)**  
- Rank interactions by **`lrscore`**  
- Select top N ligand–receptor pairs

#### Visualization
- Each dot represents a ligand–receptor interaction:
  - X-axis: ligand (source cell type)
  - Y-axis: receptor (target cell type)
  - Dot size: `lrscore` (interaction strength)
  - Color: `lr_logfc` (condition bias)

#### Interpretation

This analysis reveals:
- Which cells are sending signals
- Which cells are receiving signals
- Which interactions are:
  - Strong and consistent (`lrscore`)
  - Condition-specific (`lr_logfc`)

#### LIANA Output Metrics

| Column | Meaning | Interpretation |
|--------|--------|----------------|
| **source** | Sending cell type | Cell type producing the ligand (signal sender) |
| **target** | Receiving cell type | Cell type expressing the receptor (signal receiver) |
| **ligand_complex** | Ligand gene / complex | Signaling molecule produced by source cells |
| **receptor_complex** | Receptor gene / complex | Molecule on target cells receiving the signal |
| **lr_means** | Mean expression of ligand and receptor | Indicates whether both ligand and receptor are moderately/highly expressed |
| **expr_prod** | Product of ligand × receptor expression | Measures co-expression strength (high only if both are expressed) |
| **lr_logfc** | Log fold-change (e.g. Reg vs nonReg) | Positive = enriched in Reg, negative = enriched in nonReg |
| **cellphone_pvals** | Permutation p-value (CellPhoneDB) | Statistical significance of interaction (lower = more significant) |
| **scaled_weight** | Normalized interaction strength | Consensus-scaled interaction magnitude across methods |
| **spec_weight** | Specificity score | How unique the interaction is to specific cell–cell pairs |
| **lrscore** | LIANA consensus score (0–1) | Overall confidence of interaction across methods |
| **specificity_rank** | Rank of specificity | Lower rank = more cell-type-specific interaction |
| **magnitude_rank** | Rank of interaction strength | Lower rank = stronger interaction compared to others |


## Celltypes similarities between samples

### Method 1: Cosine similarity 

- Cosine similarity measures the **angle** between two vectors. **1** = same direction. **0** = perpendicular.

- It only cares about **pattern**, not total expression level. If all genes are 2x higher in one sample but relative pattern is the same → similarity = 1.

- It shows the **pattern** of which genes are high and low in cell type A is the same as in cell type B.

- It doesn't tell you how much higher or lower a specific gene is (e.g., 3 vs 30). Magnitude is ignored.


Cosine similarity = 1 (if all genes scaled equally) → Tells you nothing about the 10x change. Use log fold change instead.

![](figures/axt_cosine_similarity.png?v=7)

### PCA + Wasserstein Distance

##### Step 1: PCA Dimensionality Reduction

Every cell has thousands of gene expression values. PCA compresses this into just a few coordinates (PC1, PC2, etc.). Each cell becomes a single point in this simplified space, with its position representing its entire transcriptional state.

##### Step 2: Cloud Representation

For any two groups you want to compare—whether different cell types, different samples, or different conditions—you plot each group as a cloud of points in the PC space.

##### Step 3: Wasserstein Distance Calculation

The Wasserstein distance calculates how much the points in one cloud would need to move to perfectly overlap the other cloud.

- **Small distance** = clouds already overlap = groups are transcriptionally similar
- **Large distance** = clouds are far apart = groups are transcriptionally different

##### Step 4: Convert to Similarity Score

Similarity = 1 - (distance / max_distance)

- **Higher number (close to 1)** = more similar
- **Lower number (close to 0)** = more different

##### Output

A similarity score between 0 and 1 for any pair of groups, where higher means more transcriptionally similar.
 
![](figures/axt_celltype_similarity_heatmap.png?v=2)

| Aspect | Cosine Similarity | PCA + Wasserstein |
|--------|-------------------|-------------------|
| **What it reflects** | Whether the average expression pattern is the same | Whether the entire population structure is the same |
| **What it tells you** | If the average expression across all cells is similar between groups | If all cells (rare subpopulations, activated states, full heterogeneity) are similarly distributed |


### Celltypes similarities between Reg vs non Reg

#### We use subtype distribution to reveal cell type similarity between Reg vs non Reg 

###### Subcluster Composition Plot

For each subtype within a cell type:

1. Count cells from condition A
2. Count cells from condition B  
3. Calculate total = A + B
4. Proportion A = A / total
5. Proportion B = B / total
6. Red bar = Proportion A
7. Blue bar = Proportion B

Equal bars (both 0.5) = subtype equally represented in both conditions. Unequal bars = subtype enriched in one condition.

<img src="figures/axt_Fibroblast_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Osteosarcoma_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Macrophage_subcluster_composition.png?v=1" width="33%" />

<img src="figures/axt_Basophil_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Neutrophil_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Schwann_subcluster_composition.png?v=1" width="33%" />

<img src="figures/axt_T-Cells_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Chondrocyte_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Osteoblast_subcluster_composition.png?v=1" width="33%" />

<img src="figures/axt_Osteoclast_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Lymphatic_Endothelial_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_NailEpithelium_Keratinocyte_subcluster_composition.png?v=1" width="33%" />

<img src="figures/axt_SMC_Pericyte_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_B-Cells_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_OsteoProgenitor_subcluster_composition.png?v=1" width="33%" />

<img src="figures/axt_SweatGland_subcluster_composition.png?v=1" width="33%" /><img src="figures/axt_Endothelial_subcluster_composition.png?v=1" width="33%" />

###### Similarity Score Calculation

Using the proportions calculated above, we compute a similarity score for each cell type:

**Step 1:** For each subtype, calculate absolute difference = |Proportion A - Proportion B|

**Step 2:** Average all subtype differences within the cell type

**Step 3:** Similarity = 1 - (Average Difference)

**Step 4:** One bar per cell type in the similarity barplot, where bar height = similarity score

![](figures/axt_subcluster_similarity_barplot.png?v=2)

### A look into Proliferation Genes

To rule out proliferation as a driver of the Reg vs nonReg transcriptional difference, we examined expression of canonical proliferation markers.

![](figures/dotplot__AXT_Proliferation_dotplot.png?v=2)

### A look into Inflammation Genes 

To rule out inflammation as a driver of the Reg vs nonReg transcriptional difference, we examined expression of canonical inflammatory markers 

![](figures/dotplot__AXT_Inflammation_dotplot.png?v=2)

## Pathways and GO analysis 

We performed over‑representation enrichment analysis using g:Profiler similar to what we did in Osteosarcoma 

![](figures/axt_KEGG,REAC,WP_enrichment.png?v=3)

[Download here KEGG/REAC/WP enrichments](https://docs.google.com/spreadsheets/d/1u2G4ApoT3-Ru5HIvg39RoBjF9Pbp09PfZ3VQ359RIXU/edit?usp=sharing)


#### GO
![](figures/axt_GO_enrichment.png?v=3)

[Download here GO enrichments](https://docs.google.com/spreadsheets/d/1dCW49iS77fXt7F8ZQ-Q5_LRU53YpOYl0Yci_5BFxyGY/edit?usp=sharing)

#### CORUM
![](figures/axt_CORUM_enrichment.png?v=3)

[Download CORUM enrichments here](https://docs.google.com/spreadsheets/d/1a_JymRSOHq0ig11d8BZzpX7yBNqgNizO-1dcmran0to/edit?usp=sharing)

[Go to how to read go profiler results section](#how-to-read-go-profiler-results)

## MORE ANALYSIS ON THE WAY 


## REFERENCES 


---

## 📚 Software and Method Citations

This analysis pipeline is built using the following tools and methods:

- 🧬 **LIANA**  
  Dimitrov et al., 2022. *LIANA: a toolbox for cell–cell communication inference.*  
  https://github.com/saezlab/liana

- 📡 **CellPhoneDB**  
  Efremova et al., 2020. *CellPhoneDB: inferring cell–cell communication from single-cell data.*  
  https://www.cellphonedb.org/

- 🧠 **CellChat**  
  Jin et al., 2021. *CellChat: inference and analysis of cell–cell communication.*  
  https://www.nature.com/articles/s41467-021-21246-9

- 🔗 **NATMI**  
  Hou et al., 2020. *NATMI: network analysis toolkit for multilayer interactions.*  
  https://github.com/asrhou/NATMI

- 🌐 **Connectome**  
  Raredon et al., 2021. *Connectome: multi-modal analysis of cell–cell communication networks.*  
  https://www.nature.com/articles/s41598-022-07959-x

- 🧬 **SingleCellSignalR**  
  Cabello-Aguilar et al., 2020. *SingleCellSignalR: inference of intercellular communication networks.*  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7261168/

- 📊 **Scanpy**  
  Wolf et al., 2018. *Scanpy: large-scale single-cell gene expression analysis.*  
  https://scanpy.readthedocs.io/


- **gProfiler**- 
  Kolberg, L., Raudvere, U., Kuzmin, I., Adler, P., Vilo, J., & Peterson, H. (2023). g:Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Research*, 51(W1), W207–W212. https://doi.org/10.1093/nar/gkad347

