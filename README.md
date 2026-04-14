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

<img src="figures/umap_axt_Reg.png?v=1" width="600" /> <img src="figures/umap_axt_nonReg.png?v=1" width="600" />


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


## Marker genes 

## Marker Gene Annotation Strategy

To annotate cell types, I **manually selected (guessed) canonical marker genes** based on known biology from literature.

These markers were then used to validate and assign identities to clusters.

---

## Marker Gene Sets

```python
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
```
### Dotplot for marker genes 

![](figures/dotplot__axt_dotplot.png?v=1) 


### Feature plots for marker genes 

<img src="figures/umap_axt_Abcc9.png?v=1" width="33%" /><img src="figures/umap_axt_Cd79a.png?v=1" width="33%" /><img src="figures/umap_axt_Cxcr2.png?v=1" width="33%" />
<img src="figures/umap_axt_Krt14.png?v=1" width="33%" /><img src="figures/umap_axt_Krt17.png?v=1" width="33%" /><img src="figures/umap_axt_Rspo3.png?v=1" width="33%" />
<img src="figures/umap_axt_Acan.png?v=1" width="33%" /><img src="figures/umap_axt_Cd8a.png?v=1" width="33%" /><img src="figures/umap_axt_Dcn.png?v=1" width="33%" />
<img src="figures/umap_axt_Acp5.png?v=1" width="33%" /><img src="figures/umap_axt_Cdh5.png?v=1" width="33%" /><img src="figures/umap_axt_Dcstamp.png?v=1" width="33%" />
<img src="figures/umap_axt_Acta1.png?v=1" width="33%" /><img src="figures/umap_axt_Cilp2.png?v=1" width="33%" /><img src="figures/umap_axt_Emcn.png?v=1" width="33%" />
<img src="figures/umap_axt_Acta2.png?v=1" width="33%" /><img src="figures/umap_axt_Cnn1.png?v=1" width="33%" /><img src="figures/umap_axt_Eng.png?v=1" width="33%" />
<img src="figures/umap_axt_Adgre1.png?v=1" width="33%" /><img src="figures/umap_axt_Col10a1.png?v=1" width="33%" /><img src="figures/umap_axt_Erg.png?v=1" width="33%" />
<img src="figures/umap_axt_Alpl.png?v=1" width="33%" /><img src="figures/umap_axt_Col1a1.png?v=1" width="33%" /><img src="figures/umap_axt_Flt4.png?v=1" width="33%" />
<img src="figures/umap_axt_Bglap.png?v=1" width="33%" /><img src="figures/umap_axt_Col1a2.png?v=1" width="33%" /><img src="figures/umap_axt_Frzb.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd163.png?v=1" width="33%" /><img src="figures/umap_axt_Col23a1.png?v=1" width="33%" /><img src="figures/umap_axt_Gdf5.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd19.png?v=1" width="33%" /><img src="figures/umap_axt_Col2a1.png?v=1" width="33%" /><img src="figures/umap_axt_Ibsp.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd3d.png?v=1" width="33%" /><img src="figures/umap_axt_Csf1r.png?v=1" width="33%" /><img src="figures/umap_axt_Ighm.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd3e.png?v=1" width="33%" /><img src="figures/umap_axt_Cspg4.png?v=1" width="33%" /><img src="figures/umap_axt_Il2rb.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd4.png?v=1" width="33%" /><img src="figures/umap_axt_Ctsk.png?v=1" width="33%" /><img src="figures/umap_axt_Kcnj8.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd68.png?v=1" width="33%" /><img src="figures/umap_axt_Cxcl12.png?v=1" width="33%" /><img src="figures/umap_axt_Kdr.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd163.png?v=1" width="33%" /><img src="figures/umap_axt_Lepr.png?v=1" width="33%" /><img src="figures/umap_axt_Mbp.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd68.png?v=1" width="33%" /><img src="figures/umap_axt_Ly6g.png?v=1" width="33%" /><img src="figures/umap_axt_Mmp9.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd79a.png?v=1" width="33%" /><img src="figures/umap_axt_Lyve1.png?v=1" width="33%" /><img src="figures/umap_axt_Mpo.png?v=1" width="33%" />
<img src="figures/umap_axt_Cd8a.png?v=1" width="33%" /><img src="figures/umap_axt_Mrc1.png?v=1" width="33%" /><img src="figures/umap_axt_Ms4a1.png?v=1" width="33%" />
<img src="figures/umap_axt_Cdh5.png?v=1" width="33%" /><img src="figures/umap_axt_Myh11.png?v=1" width="33%" /><img src="figures/umap_axt_Pdgfra.png?v=1" width="33%" />
<img src="figures/umap_axt_Cilp2.png?v=1" width="33%" /><img src="figures/umap_axt_Pecam1.png?v=1" width="33%" /><img src="figures/umap_axt_Pdgfrb.png?v=1" width="33%" />
<img src="figures/umap_axt_Cnn1.png?v=1" width="33%" /><img src="figures/umap_axt_Pdpn.png?v=1" width="33%" /><img src="figures/umap_axt_Pi16.png?v=1" width="33%" />
<img src="figures/umap_axt_Col10a1.png?v=1" width="33%" /><img src="figures/umap_axt_Plp1.png?v=1" width="33%" /><img src="figures/umap_axt_Pmp22.png?v=1" width="33%" />
<img src="figures/umap_axt_Col1a1.png?v=1" width="33%" /><img src="figures/umap_axt_Prg4.png?v=1" width="33%" /><img src="figures/umap_axt_Prox1.png?v=1" width="33%" />
<img src="figures/umap_axt_Col1a2.png?v=1" width="33%" /><img src="figures/umap_axt_Prrx1.png?v=1" width="33%" /><img src="figures/umap_axt_Rgs5.png?v=1" width="33%" />
<img src="figures/umap_axt_Col23a1.png?v=1" width="33%" /><img src="figures/umap_axt_S100a8.png?v=1" width="33%" /><img src="figures/umap_axt_S100a9.png?v=1" width="33%" />
<img src="figures/umap_axt_Csf1r.png?v=1" width="33%" /><img src="figures/umap_axt_S100b.png?v=1" width="33%" /><img src="figures/umap_axt_Sox10.png?v=1" width="33%" />
<img src="figures/umap_axt_Sox9.png?v=1" width="33%" /><img src="figures/umap_axt_Sp7.png?v=1" width="33%" /><img src="figures/umap_axt_Tagln.png?v=1" width="33%" />
<img src="figures/umap_axt_Thy1.png?v=1" width="33%" /><img src="figures/umap_axt_Traf6.png?v=1" width="33%" /><img src="figures/umap_axt_Ucma.png?v=1" width="33%" />
<img src="figures/umap_axt_leiden.png?v=1" width="33%" /><img src="figures/umap_axt_Eng.png?v=1" width="33%" /><img src="figures/umap_axt_Erg.png?v=1" width="33%" />
<img src="figures/umap_axt_Flt4.png?v=1" width="33%" /><img src="figures/umap_axt_Kdr.png?v=1" width="33%" /><img src="figures/umap_axt_Rspo3.png?v=1" width="33%" />
<img src="figures/umap_axt_Reg.png?v=1" width="33%" />



