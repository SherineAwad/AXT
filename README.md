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

![](figures/violin_axt_preQC.png?v=1)

### Post filtering 

![](figures/violin_axt_AfterQC.png?v=1) 




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
![](figures/umap_axt.png?v=1) 

#### Per sample UMAP 

<img src="figures/umap_axt_Reg.png?v=1" width="45%" /> <img src="figures/umap_axt_nonReg?v=1" width="45%" />


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


![](figures/umap_axt_leiden.png?v=1) 


## QC per Leiden Cluster

QC metrics were visualized across Leiden clusters to assess cluster quality.
Check QC per cluster to spot and remove **low-quality or suspicious cell groups**

<img src="figures/violin_axt_QC_n_genes_by_counts.png?v=1" width="33%" /><img src="figures/violin_axt_QC_pct_counts_mt.png?v=1" width="33%" /><img src="figures/violin_axt_QC_total_counts.png?v=1" width="33%" />


