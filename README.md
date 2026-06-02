# Project Overview

This project combines our **AXT dataset** with publicly available single-cell RNA-seq samples from **GEO dataset GSE135985**.

## Our Dataset
- AXT dataset: https://github.com/SherineAwad/AXT

🔺🔺🔺 Aligned to reference genome GRCm39 and EGFP sequence is added 

## External Dataset (GSE135985)

The following samples were integrated:

- **GSM4038978** — Regen_14 DPA_scRNA-seq  
- **GSM4038982** — Uninjured_scRNA-seq  
- **GSM4038983** — Uninjured_2_scRNA-seq  
- **GSM4038986** — Non_Regen_14 DPA_scRNA-seq  

🔺🔺🔺 Raw files were pulled from GEO and re-aligned to GRCm39

## Objective

The combined dataset is used to compare regenerative versus non-regenerative conditions across injury states, enabling integrated single-cell analysis and cross-dataset comparison.

### Light filtering 

#### Pre filtering 
![](figures/violin_GSE135985_preQC.png?v=3)

#### Post filtering 
![](figures/violin_GSE135985_AfterQC.png?v=3)

### Doublet detection 

![](figures/GSE135985_scrublet_scores.png?v=3)

Estimated detectable doublet fraction = 1.6%
Overall doublet rate:
	Expected   = 6.0%
	Estimated  = 1.1%
Elapsed time: 19.4 seconds
Detected doublets: 20 / 11617
Remaining cells after filtering: 11597


### Now merge this data set with our AXT data set 

| Dataset | Cells | Genes |
|---------|-------|-------|
| Input 1 | 31318 | 27808 |
| Input 2 | 11597 | 23702 |
| Merged | 42915 | 27909 |


| Description                  | Value  |
|------------------------------|--------|
| Expected unique genes | 27909 |
| Actual genes in merge | 27909 |

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


![](figures/umap_combined_umap.png?v=6)

<img src="figures/umap_combined_nonReg.png?v=6" width="33%" /><img src="figures/umap_combined_Reg_14DPA.png?v=4" width="33%" /><img src="figures/umap_combined_Uninjured1.png?v=6" width="33%" />

<img src="figures/umap_combined_Uninjured2.png?v=6" width="33%" /><img src="figures/umap_combined_Reg.png?v=6" width="33%" /><img src="figures/umap_combined_NonReg_14DPA.png?v=6" width="33%" />

## Use run_harmony for batch effect 

##  🟢 🟢 🟢 Harmony using All Genes ~30k genes

![](figures/umap_harmony_noHVG_harmony_batch.png?v=3) 


### A look into PCAs top genes

```md

### PC0
**POS:** Tmsb4x, Srgn, Fth1, Cd74, Ctss, Plek, Fcer1g, Lcp1, Cd52, Lyz2  
**NEG:** Col1a1, Col1a2, Col3a1, Fn1, Col5a2, Sparc, Col6a3, Nedd4, Serpinh1, Cald1  

### PC1
**POS:** Igfbp7, Dcn, Bgn, Sparc, Serping1, Fstl1, Lum, Gsn, Plpp3, Igfbp5  
**NEG:** Cmss1, Lars2, Camk1d, Cdk8, Ctss, Cd74, H2-Aa, EGFP, Srgn, Fabp5  

### PC2
**POS:** Cmss1, Camk1d, Lars2, Cdk8, Gphn, Mast4, Ext1, Apbb2, Zeb1, Zbtb20  
**NEG:** Fth1, Ppia, Crip1, Ftl1, Lgals3, Lgals1, S100a6, Nme2, Gapdh, S100a11  
```

### Clustering 

![](figures/umap_combined_leiden.png?v=6)

<img src="figures/violin_combined_QC_n_genes_by_counts.png?v=6" width="33%" /><img src="figures/violin_combined_QC_total_counts.png?v=6" width="33%" /><img src="figures/violin_combined_QC_pct_counts_mt.png?v=6" width="33%" />

### Visualising marker genes 

We used the following marker genes:

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

![](figures/dotplot__combined_dotplot.png?=3)

<img src="figures/umap_combined_Nfatc1.png?v=6" width="33%" /><img src="figures/umap_combined_Col1a1.png?v=6" width="33%" /><img src="figures/umap_combined_Pecam1.png?v=6" width="33%" />

<img src="figures/umap_combined_Prrx1.png?v=6" width="33%" /><img src="figures/umap_combined_Spp1.png?v=6" width="33%" /><img src="figures/umap_combined_Cd44.png?v=6" width="33%" />

<img src="figures/umap_combined_Dcn.png?v=6" width="33%" /><img src="figures/umap_combined_Runx2.png?v=6" width="33%" /><img src="figures/umap_combined_Postn.png?v=6" width="33%" />

<img src="figures/umap_combined_Acta2.png?v=6" width="33%" /><img src="figures/umap_combined_Pdgfrb.png?v=4" width="33%" /><img src="figures/umap_combined_Ibsp.png?v=4" width="33%" />

<img src="figures/umap_combined_Pmp22.png?v=6" width="33%" /><img src="figures/umap_combined_Cspg4.png?v=4" width="33%" /><img src="figures/umap_combined_Cd79a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ctsk.png?v=6" width="33%" /><img src="figures/umap_combined_Rgs5.png?v=4" width="33%" /><img src="figures/umap_combined_Pdgfra.png?v=4" width="33%" />

<img src="figures/umap_combined_EGFP.png?v=6" width="33%" /><img src="figures/umap_combined_Cxcl12.png?v=4" width="33%" /><img src="figures/umap_combined_Mmp13.png?v=4" width="33%" />

<img src="figures/umap_combined_Col23a1.png?v=6" width="33%" /><img src="figures/umap_combined_Alpl.png?v=4" width="33%" /><img src="figures/umap_combined_Erg.png?v=4" width="33%" />

<img src="figures/umap_combined_Flt1.png?v=6" width="33%" /><img src="figures/umap_combined_Sox9.png?v=4" width="33%" /><img src="figures/umap_combined_Prg4.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt14.png?v=6" width="33%" /><img src="figures/umap_combined_Tagln.png?v=4" width="33%" /><img src="figures/umap_combined_Cd34.png?v=4" width="33%" />

<img src="figures/umap_combined_Esam.png?v=6" width="33%" /><img src="figures/umap_combined_Pdpn.png?v=4" width="33%" /><img src="figures/umap_combined_Cd68.png?v=4" width="33%" />

<img src="figures/umap_combined_Cdh5.png?v=6" width="33%" /><img src="figures/umap_combined_Emcn.png?v=4" width="33%" /><img src="figures/umap_combined_Plvap.png?v=4" width="33%" />

<img src="figures/umap_combined_Nes.png?v=6" width="33%" /><img src="figures/umap_combined_Csf1r.png?v=4" width="33%" /><img src="figures/umap_combined_Lepr.png?v=4" width="33%" />

<img src="figures/umap_combined_Rspo3.png?v=6" width="33%" /><img src="figures/umap_combined_Sp7.png?v=4" width="33%" /><img src="figures/umap_combined_Tnfrsf11a.png?v=4" width="33%" />

<img src="figures/umap_combined_Kdr.png?v=6" width="33%" /><img src="figures/umap_combined_S100a8.png?v=4" width="33%" /><img src="figures/umap_combined_Abcc9.png?v=4" width="33%" />

<img src="figures/umap_combined_S100b.png?v=6" width="33%" /><img src="figures/umap_combined_S100a9.png?v=4" width="33%" /><img src="figures/umap_combined_Mbp.png?v=4" width="33%" />

<img src="figures/umap_combined_Vwf.png?v=6" width="33%" /><img src="figures/umap_combined_Myh11.png?v=4" width="33%" /><img src="figures/umap_combined_Krt5.png?v=4" width="33%" />

<img src="figures/umap_combined_Acp5.png?v=6" width="33%" /><img src="figures/umap_combined_Cdh1.png?v=4" width="33%" /><img src="figures/umap_combined_Mrc1.png?v=4" width="33%" />

<img src="figures/umap_combined_Frzb.png?v=6" width="33%" /><img src="figures/umap_combined_Kit.png?v=4" width="33%" /> <img src="figures/umap_combined_Adgre1.png?v=4" width="33%" />

<img src="figures/umap_combined_Comp.png?v=6" width="33%" /><img src="figures/umap_combined_Acan.png?v=4" width="33%" /> <img src="figures/umap_combined_Krt17.png?v=4" width="33%" />

<img src="figures/umap_combined_Hpgds.png?v=6" width="33%" /><img src="figures/umap_combined_Kcnj8.png?v=4" width="33%" /><img src="figures/umap_combined_Dmp1.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd28.png?v=6" width="33%" /><img src="figures/umap_combined_Csf3r.png?v=4" width="33%" /><img src="figures/umap_combined_Flt4.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd4.png?v=6" width="33%" /><img src="figures/umap_combined_Hdc.png?v=4" width="33%" /><img src="figures/umap_combined_Lgr6.png?v=4" width="33%" />

<img src="figures/umap_combined_Epcam.png?v=6" width="33%" /><img src="figures/umap_combined_Pi16.png?v=4" width="33%" /><img src="figures/umap_combined_Cilp2.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd163.png?v=6" width="33%" /><img src="figures/umap_combined_Dsg3.png?v=4" width="33%" /><img src="figures/umap_combined_Col9a1.png?v=4" width="33%" />

<img src="figures/umap_combined_Plp1.png?v=6" width="33%" /><img src="figures/umap_combined_Des.png?v=4" width="33%" /><img src="figures/umap_combined_Dcstamp.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd3e.png?v=6" width="33%" /><img src="figures/umap_combined_Col2a1.png?v=4" width="33%" /><img src="figures/umap_combined_Oscar.png?v=4" width="33%" />

<img src="figures/umap_combined_Lyve1.png?v=6" width="33%" /><img src="figures/umap_combined_Ighm.png?v=4" width="33%" /><img src="figures/umap_combined_Cd3d.png?v=4" width="33%" />

<img src="figures/umap_combined_Cnn1.png?v=6" width="33%" /><img src="figures/umap_combined_Aqp5.png?v=4" width="33%" /><img src="figures/umap_combined_Scnn1a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ngfr.png?v=6" width="33%" /><img src="figures/umap_combined_Prox1.png?v=4" width="33%" /><img src="figures/umap_combined_Calcr.png?v=4" width="33%" />

<img src="figures/umap_combined_Sp8.png?v=6" width="33%" /><img src="figures/umap_combined_Sp6.png?v=4" width="33%" /><img src="figures/umap_combined_Krt19.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1g.png?v=6" width="33%" /><img src="figures/umap_combined_Ms4a1.png?v=4" width="33%" /><img src="figures/umap_combined_Pax5.png?v=4" width="33%" />

<img src="figures/umap_combined_Sox10.png?v=6" width="33%" /><img src="figures/umap_combined_Mpo.png?v=4" width="33%" /><img src="figures/umap_combined_Cd22.png?v=4" width="33%" />

<img src="figures/umap_combined_Ucma.png?v=6" width="33%" /><img src="figures/umap_combined_Cd8a.png?v=4" width="33%" /><img src="figures/umap_combined_Matn1.png?v=4" width="33%" />

<img src="figures/umap_combined_Tpsb2.png?v=6" width="33%" /><img src="figures/umap_combined_Krt7.png?v=4" width="33%" /><img src="figures/umap_combined_Ccr3.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1b.png?v=6" width="33%" /><img src="figures/umap_combined_Ms4a2.png?v=4" width="33%" /><img src="figures/umap_combined_Krt18.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt8.png?v=6" width="33%" /><img src="figures/umap_combined_Cd19.png?v=4" width="33%" /><img src="figures/umap_combined_Foxa1.png?v=4" width="33%" />

<img src="figures/umap_combined_Ly6g.png?v=6" width="33%" /><img src="figures/umap_combined_Tpsab1.png?v=4" width="33%" /><img src="figures/umap_combined_Gdf5.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd200r3.png?v=6" width="33%" /><img src="figures/umap_combined_Mcpt8.png?v=4" width="33%" /><img src="figures/umap_combined_Cpa3.png?v=4" width="33%" />

<img src="figures/umap_combined_Bglap.png?v=6" width="33%" />

### Double checking cluster doublets 

From the dotplot and feature plots, 

![](figures/umap_doublet_score.png?v=3) 



### Preliminary Annotations 

![](figures/umap_combined_celltypeON.png?v=13) 

![](figures/umap_combined_celltype.png?v=13)


# From here requires re-running -- coming soon 
### Now lets subset Fibroblast from GSE135985 and Osteosarcoma from AXT 

## Celltype Composition by Sample

| celltype      | Non_Regen_14_DPA | Reg | Regen_14_DPA | Uninjured1 | Uninjured2 | nonReg |
|---------------|------------------|-----|--------------|-------------|-------------|--------|
| Fibroblast    | 744              | 0   | 638          | 1663        | 3011        | 0      |
| Osteosarcoma  | 0                | 1108| 0            | 0           | 0           | 8940   |


![](figures/Fibroblast_Osteosarcoma_umap.png?v=2)

![](figures/Fibroblast_Osteosarcoma_perSample_umap.png?v=2)

## Celltype similarity 

### Using Pearson correlation, Spearman Correlation, and Cosine similarity 

| Method | What it asks |
|----------|-------------|
| **Pearson** | "Do the same genes increase and decrease together?" |
| **Spearman** | "Do genes keep the same rank order?" |
| **Cosine** | "Do the two groups have the same overall expression pattern?" |

### Quick interpretation

- **Pearson** → Similarity of gene expression changes.
- **Spearman** → Similarity of gene rankings.
- **Cosine** → Similarity of the overall transcriptional program.

### Main difference

- **Pearson** cares about expression values.
- **Spearman** cares about gene order/ranking.
- **Cosine** cares about expression pattern shape.

![](figures/Fibroblast_Osteosarcoma_pearson_similarity.png ?v=1)

![](figures/Fibroblast_Osteosarcoma_spearman_similarity.png?v=1)

![](figures/Fibroblast_Osteosarcoma_cosine_similarity.png?v=1)

### Using PCA (Wasserstein, MMD, and Optimal Transport) 
 
| Method | What it asks |
|----------|-------------|
| **PCA + Wasserstein** | "How different are the overall cell-state distributions between these groups?" |
| **PCA + MMD** | "Do these groups come from the same underlying cellular population?" |
| **PCA + OT (Sinkhorn/POT)** | "How easily can cells from one group be matched to cells in the other group?" |

### Quick interpretation

- **PCA + Wasserstein** → Compares where cells are located across the overall cellular landscape.
- **PCA + MMD** → Compares whether the overall population structure looks the same.
- **PCA + OT (Sinkhorn/POT)** → Compares how well cells from one group can be aligned to cells from another group.

### Main difference from Pearson / Spearman / Cosine

Pearson, Spearman and Cosine compare:

> One average expression profile vs another average expression profile.

PCA + Wasserstein, MMD and OT compare:

> Entire populations of cells vs entire populations of cells.

So they can detect differences in:
- cell-state composition
- population structure
- subpopulations
- distribution of cells across states

even when the average expression profile looks similar.


![](figures/Fibroblast_Osteosarcoma_pca_wasserstein.png?v=1)

![](figures/Fibroblast_Osteosarcoma_pca_mmd.png?v=1)

![](figures/Fibroblast_Osteosarcoma_pot_matrix.png?v=2) 

### 🚨🚨🚨 Method 7 SCOT+: SCOT+ is missing -- software issue -- contacting author 





