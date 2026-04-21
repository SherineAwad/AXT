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

---

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
    "Sweat_glands": ["Aqp5","Scnn1a","Scnn1b","Scnn1g","Krt19","Krt7","Krt8","Krt18","Krt5","Krt14","Foxa1"]

}

```

### Dotplot for marker genes 

![](figures/dotplot__axt_dotplot.png?v=6) 


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

`
## Annotations 

We used the above marker genes to annotate our celltypes as follows:

![](figures/umap_axt_celltype.png?v=5)

![](figures/umap_axt_celltypeON.png?v=5) 

🚨🚨⚠️ WARNING: Cluster 43 removed due to low-quality cells (133 cells were removed) ⚠️ 🚨🚨


## Some stats for samples and cells 

| Celltype                         | Reg   | nonReg | Total |
|----------------------------------|------:|-------:|------:|
| B-Cells                          |   259 |    291 |   550 |
| Chondrocyte                      |   185 |      3 |   188 |
| Fibroblast                       |  3903 |   2104 |  6007 |
| Lymphatic_Endothelial            |   147 |    203 |   350 |
| MSC                              |  1513 |   1076 |  2589 |
| Macrophage                       |  1774 |   3221 |  4995 |
| NailEpithelium_Keratinocyte      |   663 |    491 |  1154 |
| Neutrophil                       |   221 |    424 |   645 |
| Osteoblast                       |   250 |   1513 |  1763 |
| Osteoclast                       |   136 |    561 |   697 |
| Osteosarcoma                     |   889 |   7521 |  8410 |
| SMC_Pericyte                     |   553 |    540 |  1093 |
| Schwann                          |   176 |     54 |   230 |
| SweatGland                       |   215 |      2 |   217 |
| Synoviocyte                      |   277 |    559 |   836 |
| T-Cells                          |   523 |   1014 |  1537 |
| **Total**                        | 11684 |  19577 | 31261 |


![](figures/axt_cell_ratios.png?v=5)


###
## Differential gene expression 

### Global Regenerative vs non Regenerative 

We performed global differential gene expression analysis between Reg and non-Reg samples using the Wilcoxon rank-sum test. Genes were filtered using an adjusted p-value threshold of < 0.05. The top N genes with the highest positive and lowest negative log fold change values were selected as upregulated and downregulated genes, respectively. These genes were visualized in a heatmap displaying their global log fold change values (Reg vs non-Reg). The heatmap represents differential expression effect sizes and is used for visualization only, not for statistical testing.

### N=20

![](figures/axt_heatmap.png?v=5)

#### List of DGE per sample with adj-pvalue \< 0.05 is in the link below

[Download here DGE per sample `<0.05`](https://docs.google.com/spreadsheets/d/1hDyW7gU68zLsds5Pb64CONviXox_iGGnQTMFCjDpznE/edit?usp=sharing)

### Now deeper look into celltypes

We performed differential gene expression analysis stratified by cell type to compare Reg and nonReg conditions. For each cell type, cells were subsetted and a Wilcoxon rank-sum test was applied using Scanpy’s `rank_genes_groups` function with `sample` as the grouping variable and nonReg set as the reference. This approach identifies genes that are differentially expressed in Reg relative to nonReg within each individual cell type, rather than across the full dataset. For each cell type, we extracted log fold-changes and adjusted p-values, filtered for statistically significant genes (adjusted p-value < threshold), and removed duplicate gene entries by retaining the strongest signal. The top N upregulated and downregulated genes per cell type were selected based on log fold-change and visualised in a heatmap summarising Reg versus nonReg effects across all cell types.

### N=5 per celltype 

![](figures/axt_celltype_heatmap.png?v=6)

#### List of DGE with adj-pvalue \< 0.05 is in the link below

[Download here DGE per celltype `<0.05`](https://docs.google.com/spreadsheets/d/1cb3aBfwGRP9NWf71YrZzLOhstrRxBVJ3vmAqKtGivU4/edit?usp=sharing)

###### ⚠️ ⚠️ I have the full list of DGE, but the file is too large to upload here


## 🧬 Cell–Cell Interaction Analysis (LIANA)

This pipeline identifies potential **ligand–receptor interactions** from single-cell RNA-seq data using the LIANA framework. LIANA integrates multiple established methods into a single consensus scoring system.

#### Top lrScores interactions with minimum |logFC| 0.5

![](figures/axt_liana_dotplot.png?=3)


#### Click the link below for full list of cell-cell interactions

[Download here cell to cell interactions](https://docs.google.com/spreadsheets/d/1CnfeNbDfsNKgrBq6TPma7-F5rTUY8mQ2vthxwdHB3ng/edit?usp=sharing)

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

![](figures/axt_cosine_similarity.png?v=3) 


## Pathways and GO analysis 

We performed over‑representation enrichment analysis using g:Profiler against the following functional annotation databases:

- **KEGG**: Manually curated pathway maps for metabolism, signaling, and diseases.
- **REAC (Reactome)**: Peer‑reviewed, hierarchical pathway database covering molecular processes.
- **GO (Gene Ontology)**: Three structured ontologies (BP, CC, MF) describing gene function.
- **WP (WikiPathways)**: Community‑curated, open platform for biological pathways.
- **TFs (TRANSFAC)**: Transcription factor binding sites and regulated targets.
- **CORUM**: Experimentally validated mammalian protein complexes.

For each cell type, we split differentially expressed genes into UP‑regulated (logFC > 0) and DOWN‑regulated (logFC < 0) lists and tested each against the above databases.  
Significance was assessed with a one‑sided Fisher’s exact test, and p‑values were corrected for multiple testing using the g:SCS method with an adjusted cutoff of **0.01**.  
Within each cell type and direction, we retained the top **10** most significant terms (lowest raw p‑value).  

In the dotplot, the **number of overlapping genes** between the query list and a term is represented by the **dot shape**:  
- Circle → 1–5 genes  
- Square → 6–10 genes  
- Triangle → 11–15 genes  
- Diamond → ≥16 genes  

Dot **colour** indicates the significance level (–log₁₀ adjusted p‑value). The plot is split into two panels: UP‑regulated (left) and DOWN‑regulated (right).


#### KEGG
![](figures/axt_enrichment_KEGG_dotplot.png?v=2)


#### REAC
![](figures/axt_enrichment_REAC_dotplot.png?v=2)


#### WP
![](figures/axt_enrichment_WP_dotplot.png?v=2)


##### GO
![](figures/axt_enrichment_GO_dotplot.png?v=2)


##### TFs 
![](figures/axt_enrichment_TF_dotplot.png?v=2)


#### CORUM
![](figures/axt_enrichment_CORUM_dotplot.png?v=2)


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

