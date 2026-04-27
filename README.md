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

## 🚨⚠️ **The deconatimantion step will affect some of the downstream analysis, so please hold on till I update to reflect the above and remove this sign** ⚠️  🚨
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

![](figures/ Osteosarcoma_volcano.png?v=6) 

#### Dotplot for genes of interest in Osteosarcoma 

![](figures/dotplot__Osteosarcoma_dotplot.png?v=1)

#### Violin plots for genes of interest in Osteosarcoma 

<img src="figures/violin_Osteosarcoma_Aspn_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Bglap2_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Bglap_violin.png?v=1" width="33%" />
<img src="figures/violin_Osteosarcoma_Kazald1_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Nbl1_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Notum_violin.png?v=1" width="33%" />
<img src="figures/violin_Osteosarcoma_Panx3_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Alpl_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Cdkn1a_violin.png?v=1" width="33%" />
<img src="figures/violin_Osteosarcoma_Dcn_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Omd_violin.png?v=1" width="33%" /><img src="figures/violin_Osteosarcoma_Sparc_violin.png?v=1" width="33%" />
<img src="figures/violin_Osteosarcoma_Wif1_violin.png?v=1" width="33%" />


#### Feature plots for genes of interest in Osteosarcoma

<img src="figures/Osteosarcoma_Notum_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Bglap2_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Aspn_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Bglap_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Nbl1_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Panx3_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Kazald1_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Dcn_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Wif1_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Cdkn1a_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Alpl_featureplot.png?v=1" width="33%" /><img src="figures/Osteosarcoma_Omd_featureplot.png?v=1" width="33%" />
<img src="figures/Osteosarcoma_Sparc_featureplot.png?v=1" width="33%" />


## Pathways and GO Enrichment for  Osteosarcoma 

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
![](figures/Osteosarcoma_KEGG,REAC,WP_enrichment.png?v=1)

#### Using GO
![](figures/Osteosarcoma_GO_enrichment.png?v-1) 

#### Using CORUM 
![](figures/Osteosarcoma_CORUM_enrichment.png?v=1)  


## Zooming on Fibroblast 

We subset Fibroblast cell type, below is a summary:

| Metric   | Value |
|----------|------:|
| Cells    | 3,586 |
| Genes    | 27,808 |

| Sample  | Count |
|---------|------:|
| Reg     | 2,759 |
| nonReg  |   827 |

![](figures/Fibroblast_umap.png?v=1)

![](figures/Fibroblast_perSample_umap.png?v=1)

### Differential gene expression in Fibroblast

Similar to what we did with Osteosarcoma, Top N genes are selected by first sorting all genes by adjusted p-value in ascending order (most statistically significant first), with log fold-change used only as a secondary criterion to break ties while preserving directionality. The top-ranked genes are then taken globally from this ordered list.

#### Top 100 
![](figures/Fibroblast_heatmap.png?v=3)

Now overlapping with Secreted genes 
  
#### Top 10 genes only are labeled 
![](figures/Fibroblast_volcano.png?v=4)

[Click here to download Fibroblast DGE with adj-pvalue \< 0.05 ](https://docs.google.com/spreadsheets/d/1h8QIfU7iIcgeYpKBaU6IiNVt5lNo3TAMVOqjLGB6V2k/edit?usp=sharing)

#### Dotplot for genes of interest in Fibroblast  

![](figures/Fibroblast_dotplot.png?v=2)

#### Violin plots for genes of interest in Fibroblast 

![](figures/Fibroblast_violin_Bmp2.png?v=3) 

![](figures/Fibroblast_violin_Bmp4.png?v=3)

#### Feature plots for genes of interest in Fibroblast 


<img src="figures/Fibroblast_feature_Bmp2_Reg.png?v=1" width="45%" /> <img src="figures/Fibroblast_feature_Bmp2_nonReg.png?v=1" width="45%" />
<img src="figures/Fibroblast_feature_Bmp4_Reg.png?v=1" width="45%" /> <img src="figures/Fibroblast_feature_Bmp4_nonReg.png?v=1" width="45%" />


## Patways and GO enrichments for Fibroblast  

We followed the same approach as in Osteosarcoma 

#### Using KEGG, REACTOME, and WP
![](figures/Fibroblast_KEGG,REAC,WP_enrichment.png?v=1)

#### Using GO
![](figures/Fibroblast_GO_enrichment.png?v-1) 

#### Using CORUM
![](figures/Fibroblast_CORUM_enrichment.png?v=1)   


## Now zooming on OsteoProgenitor 

We subset the OsteoProgenitor and a summary is below: 

| Metric   | Value |
|----------|------:|
| Cells    | 3,138 |
| Genes    | 27,808 |

| Sample  | Count |
|---------|------:|
| nonReg  | 1,797 |
| Reg     | 1,341 |

![](figures/OsteoProgenitor_umap.png?v=1)

![](figures/OsteoProgenitor_perSample_umap.png?v=1)


### Differential gene expression in OsteoProgenitor 

We followed the same approach we did with Osteosarcoma (same script)

#### Top 100 genes
![](figures/OsteoProgenitor_heatmap.png?v=4)

#### Top 10 genes only are labeled
![](figures/OsteoProgenitor_volcano.png?v=4)

#### List of DGE with adj-pvalue \< 0.05 is in the link below
[Click here to download DGE with adj-pvalue \< 0.05 in OsteoProgenitor](https://docs.google.com/spreadsheets/d/1VDBnjvt99mxZ49HypUTrHbFizVHImcwiMTfwJaJ4dUg/edit?usp=sharing)

### Global Regenerative vs non Regenerative
 
Here, we followed the same approach we did in Osteosarcoma for dge, but we compare Reg vs non Reg. 

#### Top 100 
![](figures/axt_heatmap.png?v=4)

#### Top 10
![](figures/axt_volcano.png?v=3)

#### List of DGE with adj-pvalue \< 0.05 is in the link below
[Click here to download DGE with adj-pvalue \< 0.05 in Reg vs non Reg](https://docs.google.com/spreadsheets/d/16wIoiAhR0IzaHF47odNBZvvNrYW2A_517pXosfPgD70/edit?usp=sharing)

 

### Now look into celltypes

We performed differential gene expression (DGE) analysis to identify genes differentially expressed across cell types. For each cell type, we applied a Wilcoxon rank-sum test using Scanpy’s rank_genes_groups function with celltype as the grouping variable, comparing each cell type against all remaining cells (one-vs-rest approach). We extracted log fold-changes and adjusted p-values for each comparison, and filtered results using an adjusted p-value threshold of 0.05 to retain statistically significant genes. For downstream visualisation, we selected the top 10 upregulated and top 10 downregulated genes per cell type (--N 10), ranking genes primarily by adjusted p-value with log fold-change used as a secondary criterion within each direction. The final results were visualised in a heatmap displaying log fold-changes of the selected genes across all cell types.

### N=10 per celltype 

![](figures/axt_celltype_heatmap.png?v=9)

#### List of DGE with adj-pvalue \< 0.05 is in the link below

[Download here DGE per celltype `<0.05`](https://docs.google.com/spreadsheets/d/1dD6gh6XWCSjBGMEb5zPAmC-FElPK37NdPfCnGyj8vqE/edit?usp=sharing)

### 🚨🚨🚨 I have the full list of DGE, but the file is too large to upload here


## 🧬 Cell–Cell Interaction Analysis (LIANA)

This pipeline identifies potential **ligand–receptor interactions** from single-cell RNA-seq data using the LIANA framework. LIANA integrates multiple established methods into a single consensus scoring system.

#### Top lrScores interactions with minimum |logFC| 0.0 (Mostly rank by lrScore)

![](figures/axt_liana_dotplot.png?=6)


#### Click the link below for full list of cell-cell interactions

[Download here cell to cell interactions](https://drive.google.com/file/d/16D9fNklX820b8p9aXMq-Yr8vVWbLPOqV/view?usp=sharing) 

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

![](figures/axt_cosine_similarity.png?v=6) 


## Pathways and GO analysis 

We performed over‑representation enrichment analysis using g:Profiler similar to what we did in Osteosarcoma 

#### GO
![](figures/axt_GO_enrichment.png?v=1)

#### KEGG, REAC, WP

![](figures/axt_KEGG,REAC,WP_enrichment.png?v=1)

#### CORUM
![](figures/axt_CORUM_enrichment.png?v=1)


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

