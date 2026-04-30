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


### Clustering 

![](figures/umap_combined_leiden.png?v=1)

<img src="figures/violin_combined_QC_n_genes_by_counts.png?v=1" width="33%" /><img src="figures/violin_combined_QC_total_counts.png?v=1" width="33%" /><img src="figures/violin_combined_QC_pct_counts_mt.png?v=1" width="33%" />


### Visualising marker genes 


![](figures/dotplot__combined_dotplot.png?=1)


<img src="figures/umap_combined_Col1a1.png?v=1" width="33%" /><img src="figures/umap_combined_Prrx1.png?v=1" width="33%" /><img src="figures/umap_combined_Dcn.png?v=1" width="33%" />
<img src="figures/umap_combined_Cd44.png?v=1" width="33%" /><img src="figures/umap_combined_Spp1.png?v=1" width="33%" /><img src="figures/umap_combined_Nfatc1.png?v=1" width="33%" />
<img src="figures/umap_combined_Rgs5.png?v=1" width="33%" /><img src="figures/umap_combined_Runx2.png?v=1" width="33%" /><img src="figures/umap_combined_Pecam1.png?v=1" width="33%" />
<img src="figures/umap_combined_Postn.png?v=1" width="33%" /><img src="figures/umap_combined_Acta2.png?v=1" width="33%" /><img src="figures/umap_combined_Pmp22.png?v=1" width="33%" />
<img src="figures/umap_combined_Ibsp.png?v=1" width="33%" /><img src="figures/umap_combined_Ctsk.png?v=1" width="33%" /><img src="figures/umap_combined_leiden.png?v=1" width="33%" />
<img src="figures/umap_combined_Mmp13.png?v=1" width="33%" /><img src="figures/umap_combined_Cspg4.png?v=1" width="33%" /><img src="figures/umap_combined_Cxcl12.png?v=1" width="33%" />
<img src="figures/umap_combined_Pdgfrb.png?v=1" width="33%" /><img src="figures/umap_combined_Sox9.png?v=1" width="33%" /><img src="figures/umap_combined_Prg4.png?v=1" width="33%" />
<img src="figures/umap_combined_Flt1.png?v=1" width="33%" /><img src="figures/umap_combined_Krt14.png?v=1" width="33%" /><img src="figures/umap_combined_Tagln.png?v=1" width="33%" />
<img src="figures/umap_combined_Pdgfra.png?v=1" width="33%" /><img src="figures/umap_combined_Erg.png?v=1" width="33%" /><img src="figures/umap_combined_Alpl.png?v=1" width="33%" />
<img src="figures/umap_combined_Pdpn.png?v=1" width="33%" /><img src="figures/umap_combined_Col23a1.png?v=1" width="33%" /><img src="figures/umap_combined_Plvap.png?v=1" width="33%" />
<img src="figures/umap_combined_Cd34.png?v=1" width="33%" /><img src="figures/umap_combined_Esam.png?v=1" width="33%" /><img src="figures/umap_combined_Cd68.png?v=1" width="33%" />
<img src="figures/umap_combined_Nes.png?v=1" width="33%" /><img src="figures/umap_combined_Cdh5.png?v=1" width="33%" /><img src="figures/umap_combined_Csf1r.png?v=1" width="33%" />
<img src="figures/umap_combined_Sp7.png?v=1" width="33%" /><img src="figures/umap_combined_Lepr.png?v=1" width="33%" /><img src="figures/umap_combined_Emcn.png?v=1" width="33%" />
<img src="figures/umap_combined_Rspo3.png?v=1" width="33%" /><img src="figures/umap_combined_S100a9.png?v=1" width="33%" /><img src="figures/umap_combined_S100a8.png?v=1" width="33%" />
<img src="figures/umap_combined_Acp5.png?v=1" width="33%" /><img src="figures/umap_combined_Mbp.png?v=1" width="33%" /><img src="figures/umap_combined_Myh11.png?v=1" width="33%" />
<img src="figures/umap_combined_Kdr.png?v=1" width="33%" /><img src="figures/umap_combined_Krt5.png?v=1" width="33%" /><img src="figures/umap_combined_Tnfrsf11a.png?v=1" width="33%" />
<img src="figures/umap_combined_Mrc1.png?v=1" width="33%" /><img src="figures/umap_combined_Acan.png?v=1" width="33%" /><img src="figures/umap_combined_Bglap.png?v=1" width="33%" />
<img src="figures/umap_combined_Comp.png?v=1" width="33%" /><img src="figures/umap_combined_Abcc9.png?v=1" width="33%" /><img src="figures/umap_combined_Dmp1.png?v=1" width="33%" />
<img src="figures/umap_combined_Frzb.png?v=1" width="33%" /><img src="figures/umap_combined_S100b.png?v=1" width="33%" /><img src="figures/umap_combined_Vwf.png?v=1" width="33%" />
<img src="figures/umap_combined_Krt17.png?v=1" width="33%" /><img src="figures/umap_combined_Hdc.png?v=1" width="33%" /><img src="figures/umap_combined_Kit.png?v=1" width="33%" />
<img src="figures/umap_combined_Pi16.png?v=1" width="33%" /><img src="figures/umap_combined_Cdh1.png?v=1" width="33%" /><img src="figures/umap_combined_Cd79a.png?v=1" width="33%" />
<img src="figures/umap_combined_Epcam.png?v=1" width="33%" /><img src="figures/umap_combined_Lgr6.png?v=1" width="33%" /><img src="figures/umap_combined_Adgre1.png?v=1" width="33%" />
<img src="figures/umap_combined_Kcnj8.png?v=1" width="33%" /><img src="figures/umap_combined_Col2a1.png?v=1" width="33%" /><img src="figures/umap_combined_Hpgds.png?v=1" width="33%" />
<img src="figures/umap_combined_Cnn1.png?v=1" width="33%" /><img src="figures/umap_combined_Lyve1.png?v=1" width="33%" /><img src="figures/umap_combined_Cd4.png?v=1" width="33%" />
<img src="figures/umap_combined_Csf3r.png?v=1" width="33%" /><img src="figures/umap_combined_Flt4.png?v=1" width="33%" /><img src="figures/umap_combined_Cd28.png?v=1" width="33%" />
<img src="figures/umap_combined_Cilp2.png?v=1" width="33%" /><img src="figures/umap_combined_Prox1.png?v=1" width="33%" /><img src="figures/umap_combined_Cd163.png?v=1" width="33%" />
<img src="figures/umap_combined_Col9a1.png?v=1" width="33%" /><img src="figures/umap_combined_Scnn1a.png?v=1" width="33%" /><img src="figures/umap_combined_Des.png?v=1" width="33%" />
<img src="figures/umap_combined_Dcstamp.png?v=1" width="33%" /><img src="figures/umap_combined_Dsg3.png?v=1" width="33%" /><img src="figures/umap_combined_Aqp5.png?v=1" width="33%" />
<img src="figures/umap_combined_Cd3e.png?v=1" width="33%" /><img src="figures/umap_combined_Plp1.png?v=1" width="33%" /><img src="figures/umap_combined_Sp6.png?v=1" width="33%" />
<img src="figures/umap_combined_Ngfr.png?v=1" width="33%" /><img src="figures/umap_combined_Cd3d.png?v=1" width="33%" /><img src="figures/umap_combined_Krt18.png?v=1" width="33%" />
<img src="figures/umap_combined_Krt8.png?v=1" width="33%" /><img src="figures/umap_combined_Cd200r3.png?v=1" width="33%" /><img src="figures/umap_combined_Sp8.png?v=1" width="33%" />
<img src="figures/umap_combined_Ms4a1.png?v=1" width="33%" /><img src="figures/umap_combined_Cd22.png?v=1" width="33%" /><img src="figures/umap_combined_Sox10.png?v=1" width="33%" />
<img src="figures/umap_combined_Oscar.png?v=1" width="33%" /><img src="figures/umap_combined_Scnn1b.png?v=1" width="33%" /><img src="figures/umap_combined_Krt19.png?v=1" width="33%" />
<img src="figures/umap_combined_Cd8a.png?v=1" width="33%" /><img src="figures/umap_combined_Ccr3.png?v=1" width="33%" /><img src="figures/umap_combined_Krt7.png?v=1" width="33%" />
<img src="figures/umap_combined_Ucma.png?v=1" width="33%" /><img src="figures/umap_combined_Pax5.png?v=1" width="33%" /><img src="figures/umap_combined_Mpo.png?v=1" width="33%" />
<img src="figures/umap_combined_Cd19.png?v=1" width="33%" /><img src="figures/umap_combined_Mcpt8.png?v=1" width="33%" /><img src="figures/umap_combined_Scnn1g.png?v=1" width="33%" />
<img src="figures/umap_combined_Cpa3.png?v=1" width="33%" /><img src="figures/umap_combined_Gdf5.png?v=1" width="33%" /><img src="figures/umap_combined_Ms4a2.png?v=1" width="33%" />
<img src="figures/umap_combined_Ly6g.png?v=1" width="33%" />


