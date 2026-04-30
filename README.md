# Project Overview

This project combines our **AXT dataset** with publicly available single-cell RNA-seq samples from **GEO dataset GSE135985**.

## Our Dataset
- AXT dataset: https://github.com/SherineAwad/AXT

## External Dataset (GSE135985)

The following samples were integrated:

- **GSM4038978** — Regen_14 DPA_scRNA-seq  
- **GSM4038982** — Uninjured_scRNA-seq  
- **GSM4038983** — Uninjured_2_scRNA-seq  
- **GSM4038986** — Non_Regen_14 DPA_scRNA-seq  

## Objective

The combined dataset is used to compare regenerative versus non-regenerative conditions across injury states, enabling integrated single-cell analysis and cross-dataset comparison.

### Light filtering 

#### Pre filtering 
![](figures/violin_GSE135985_preQC.png?v=1)

#### Post filtering 
![](figures/violin_GSE135985_AfterQC.png?v=1)

### Doublet detection 

![](figures/GSE135985_scrublet_scores.png?v=1)

Estimated detectable doublet fraction = 5.6%
Overall doublet rate:
	Expected   = 6.0%
	Estimated  = 0.4%
Elapsed time: 26.1 seconds
Detected 7 doublets (0.1%)

7 doublets were detected and removed. 


### Now merge this data set with our AXT data set 

| Dataset      | Cells  | Genes  |
|--------------|--------|--------|
| AXT          | 31,318 | 27,808 |
| GSE135985    | 9,986  | 18,584 |
| Shared genes | —      | 17,084 |
| Merged       | 41,304 | 17,084 |


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


![](figures/umap_combined_umap.png?v=1)

<img src="figures/umap_combined_nonReg.png?v=1" width="33%" /><img src="figures/umap_combined_Reg_14DPA.png?v=1" width="33%" /><img src="figures/umap_combined_Uninjured1.png?v=1" width="33%" />

<img src="figures/umap_combined_Uninjured2.png?v=1" width="33%" /><img src="figures/umap_combined_Reg.png?v=1" width="33%" /><img src="figures/umap_combined_NonReg_14DPA.png?v=1" width="33%" />

### Now batch correction using harmony 

![](figures/umap_combined_harmony_sample.png?v=1) 

Biologically distinct samples !!



