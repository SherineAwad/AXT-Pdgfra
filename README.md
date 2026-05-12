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
![](figures/violin_GSE135985_preQC.png?)

#### Post filtering 
![](figures/violin_GSE135985_AfterQC.png?)

### Doublet detection 

![](figures/GSE135985_scrublet_scores.png?)

Estimated detectable doublet fraction = 5.6%
Overall doublet rate:
	Expected   = 6.0%
	Estimated  = 0.4%
Elapsed time: 26.1 seconds
Detected 7 doublets (0.1%)

7 doublets were detected and removed. 


### Now merge this data set with our AXT data set 
| Dataset   | Cells | Genes |
|--------------|--------|--------|
| AXT          | 31,318 | 27,808 |
| GSE135985    | 9,986  | 18,584 |
| Merged       | 41,304 | 29,308 |

| Description                  | Value  |
|------------------------------|--------|
| Expected unique genes        | 29,308 |
| Actual genes in merge        | 29,308 |


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


![](figures/umap_combined_umap.png?v=3)

<img src="figures/umap_combined_nonReg.png?v=3" width="33%" /><img src="figures/umap_combined_Reg_14DPA.png?v=3" width="33%" /><img src="figures/umap_combined_Uninjured1.png?v=3" width="33%" />

<img src="figures/umap_combined_Uninjured2.png?v=3" width="33%" /><img src="figures/umap_combined_Reg.png?v=3" width="33%" /><img src="figures/umap_combined_NonReg_14DPA.png?v=3" width="33%" />

## Use run_harmony for batch effect 

##  🟢 🟢 🟢 Harmony using All Genes ~30k genes

![](figures/umap_harmony_noHVG_harmony_batch.png?v=1) 


![](figures/PCA_harmony_noHVG_pc1_pc2_batch.png?v=1)

### PC1
- Positive genes: Cmss1, Camk1d, Lars2, Rack1, Col1a1, Gphn, Cdk8, Fn1, Col1a2, EGFP  
- Negative genes: Tmsb4x, Fth1, Ftl1, H2afz, Gnb2l1, Cst3, Crip1, Shfm1, Igfbp7, Gm8730  
---

### PC2
- Positive genes: Sparc, Col3a1, Igfbp7, Bgn, Col1a2, Dcn, Serpinh1, Col1a1, Fstl1, Lum  
- Negative genes: Cmss1, Srgn, Camk1d, Cdk8, Plek, Cd74, Ctss, Lcp1, Ptprc, Lars2  
---

### PC3
- Positive genes: Igfbp7, Dcn, Serping1, Ebf1, Sparc, Plpp3, Igfbp5, Zbtb20, Gsn, Bgn  
- Negative genes: Fabp5, Pfn1, H2az1, Ctss, Mif, Rbm3, Ppia, Actg1, Cd74, Slc25a5  

---

## 🟢 🟢 🟢 Harmony using HVGs ~2000 genes

![](figures/umap_harmony_harmony_batch.png?v=1) 

![](figures/PCA_harmony_pc1_pc2_batch.png?v=1) 

### PC1
- Positive genes: Ctss, Cd74, Fcer1g, Plek, Lyz2, Cd52, Tyrobp, H2-Ab1, Ptprc, Il1b  
- Negative genes: Dcn, Igfbp5, Serpine2, Apod, Aspn, Crispld2, Sparcl1, Cxcl12, Rgs5, Aqp1  
---

### PC2
- Positive genes: Dcn, Cst3, Apod, Igfbp5, Ccl2, Apoe, Ccl7, Ctsk, Cxcl12, Aspn  
- Negative genes: Kcnq5, Taco1, Ibsp, Epha3, Phlpp1, Spp1, Lgals7, Fabp5, Sox6, Ly6d  
---

### PC3
- Positive genes: Lgals7, Rgs5, Apoe, Sparcl1, Krtdap, Krt14, Cst3, Ly6d, Perp, Gm42418  
- Negative genes: Spp1, Tnc, Kcnq5, Ibsp, Taco1, Glis3, Epha3, Phlpp1, Col11a1, Gm20186  


## Explanation: Separation (All Genes) vs Overlap (HVGs)

### 1. Full-gene PCA → Separation

#### What drives the structure:
- **EGFP** (batch-specific transgene)
  - Present only in one batch, appears in PC1 positive genes
- Strong **ECM / fibroblast programs**
  - Col1a1, Col1a2, Fn1, Sparc, Dcn, Bgn
- **Global RNA / housekeeping variation**
  - Fth1, Ftl1, H2afz
- Broad metabolic / structural shifts

#### What this means:
- PCA is dominated by **global transcriptional state + batch-specific markers**
- EGFP contributes to PC1 alongside ECM and structural genes
- Signals are distributed across many low-variance genes (including EGFP, which has low within-batch variance but high between-batch difference)
- Strong tissue-level and batch-specific differences emerge

#### Interpretation:
> Separation between datasets reflects differences in **tissue architecture, composition, global expression state, and the presence of batch-specific markers like EGFP**.
---

### 2. HVG PCA → Overlap

#### What drives the structure:
- **Immune programs**
  - Cd74, Ptprc, Lyz2, Tyrobp, Il1b
- **Stromal programs**
  - Dcn, Igfbp5, Cxcl12, Apoe
- **Epithelial / lineage programs**
  - Krt14, Spp1, Rgs5

#### What this means:
- Only high-variance, identity-defining genes are kept
- **EGFP is NOT among the HVGs** — it has low within-batch variance (all cells in batch A express it uniformly; all cells in batch B lack it), so it is excluded by HVG selection
- Global low-variance structure (ECM tone, housekeeping, EGFP) is removed
- PCA reflects **cell-type relationships only**

#### Interpretation:
> Overlap appears because datasets share the same **cell identity structure**, once global state variation and batch-specific markers like EGFP are removed. EGFP's absence from HVGs explains why the separation disappears.
---

### 📈 Correlation between PCA and EGFP

##### Pearson Correlation (r)

r = cov(EGFP, PC1) / (σ_EGFP × σ_PC1)

**What it calculates:**

For each cell: Is EGFP above its mean WHEN PC1 is above its mean?

- Same direction → positive contribution
- Opposite direction → negative contribution
- Sum across all cells → divide by spread of both variables

**Result: -1 to +1**

##### Correlation Results

**EGFP-PC1 correlation: 0.64**

Strong positive relationship. EGFP drives PC1. Cells with EGFP sit at one end of PC1. Cells without EGFP sit at the other end.

**EGFP-PC2 correlation: 0.12**

No meaningful relationship. EGFP does not drive PC2. PC2 is driven by other biology.

**Conclusion:**

EGFP is a contributor to batch separation along PC1 in full-gene PCA (r = 0.64). Other genes (ECM, etc) also contribute. Remove EGFP via HVG selection → separation disappears, revealing shared cell identity.


### 📈 Quantifying EGFP contribution to PCA separation

To assess how much EGFP contributes to the separation between batches, we performed a linear regression of **PC1 and PC2 scores** on EGFP expression. The results are:

| Principal Component | R² (variance explained by EGFP) |
|-------------------|--------------------------------|
| PC1               | 0.404                          |
| PC2               | 0.014                          |

**Interpretation:**

- **PC1:** EGFP explains ~40% of the variance along PC1. This indicates that EGFP is a **significant contributor** to the separation observed in PC1, but other genes also contribute.  
- **PC2:** EGFP explains only ~1.4% of the variance along PC2. This shows that EGFP has **minimal influence** on PC2, which is likely driven by other transcriptional programs or cell-type-specific variation.

**Conclusion:** EGFP is **one of several contributors** to batch separation in the PCA of all genes, consistent with the observation that highly variable genes (HVG) PCA — which excludes EGFP — shows overlapping batches.

![](figures/combined_noHVG_egfp_pc_r2.png?v=1) 


### 🧬🧬 Key Insight

#### Full genes:
- Capture **global biological state + batch-specific markers (EGFP) + tissue composition**
- EGFP appears in PC1 positive genes → contributes to separation
- → leads to strong separation

#### HVGs:
- **EGFP is excluded** (low within-batch variance despite high expression)
- Capture **cell identity structure only**
- → leads to overlap
---

### Final Interpretation

> The datasets are biologically consistent at the level of cell types, but differ in global transcriptional state, tissue composition, and the presence of a batch-specific marker (EGFP). 
> 
> **Evidence for EGFP as a contributor:**
> - EGFP appears in PC1 positive genes in full-gene PCA
> - EGFP is **absent from the HVG list** (low within-batch variance)
> - Separation disappears when EGFP and other low-variance global features are removed (HVG analysis)
> 
> HVG analysis deliberately excludes EGFP and other low-variance genes, revealing shared cellular identity that is otherwise masked by batch-specific and global state differences.



### Why EGFP is not in HVG 

![](figures/combined_noHVG_EGFP_hist.png?v=1)

![](figures/combined_noHVG_EGFP_violin.png?v=1) 

![](figures/combined_noHVG_EGFP_umap.png?v=1) 



### WILL COME BACK TO THIS 


### Clustering 

![](figures/umap_combined_leiden.png?v=4)

<img src="figures/violin_combined_QC_n_genes_by_counts.png?v=4" width="33%" /><img src="figures/violin_combined_QC_total_counts.png?v=4" width="33%" /><img src="figures/violin_combined_QC_pct_counts_mt.png?v=4" width="33%" />

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


![](figures/dotplot__combined_dotplot.png?=2)

<img src="figures/umap_combined_Nfatc1.png?v=4" width="33%" /><img src="figures/umap_combined_Col1a1.png?v=4" width="33%" /><img src="figures/umap_combined_Pecam1.png?v=4" width="33%" />

<img src="figures/umap_combined_Prrx1.png?v=4" width="33%" /><img src="figures/umap_combined_Spp1.png?v=4" width="33%" /><img src="figures/umap_combined_Cd44.png?v=4" width="33%" />

<img src="figures/umap_combined_Dcn.png?v=4" width="33%" /><img src="figures/umap_combined_Runx2.png?v=4" width="33%" /><img src="figures/umap_combined_Postn.png?v=4" width="33%" />

<img src="figures/umap_combined_Acta2.png?v=4" width="33%" /><img src="figures/umap_combined_Pdgfrb.png?v=4" width="33%" /><img src="figures/umap_combined_Ibsp.png?v=4" width="33%" />

<img src="figures/umap_combined_Pmp22.png?v=4" width="33%" /><img src="figures/umap_combined_Cspg4.png?v=4" width="33%" /><img src="figures/umap_combined_Cd79a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ctsk.png?v=4" width="33%" /><img src="figures/umap_combined_Rgs5.png?v=4" width="33%" /><img src="figures/umap_combined_Pdgfra.png?v=4" width="33%" />

<img src="figures/umap_combined_EGFP.png?v=4" width="33%" /><img src="figures/umap_combined_Cxcl12.png?v=4" width="33%" /><img src="figures/umap_combined_Mmp13.png?v=4" width="33%" />

<img src="figures/umap_combined_Col23a1.png?v=4" width="33%" /><img src="figures/umap_combined_Alpl.png?v=4" width="33%" /><img src="figures/umap_combined_Erg.png?v=4" width="33%" />

<img src="figures/umap_combined_Flt1.png?v=4" width="33%" /><img src="figures/umap_combined_Sox9.png?v=4" width="33%" /><img src="figures/umap_combined_Prg4.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt14.png?v=4" width="33%" /><img src="figures/umap_combined_Tagln.png?v=4" width="33%" /><img src="figures/umap_combined_Cd34.png?v=4" width="33%" />

<img src="figures/umap_combined_Esam.png?v=4" width="33%" /><img src="figures/umap_combined_Pdpn.png?v=4" width="33%" /><img src="figures/umap_combined_Cd68.png?v=4" width="33%" />

<img src="figures/umap_combined_Cdh5.png?v=4" width="33%" /><img src="figures/umap_combined_Emcn.png?v=4" width="33%" /><img src="figures/umap_combined_Plvap.png?v=4" width="33%" />

<img src="figures/umap_combined_Nes.png?v=4" width="33%" /><img src="figures/umap_combined_Csf1r.png?v=4" width="33%" /><img src="figures/umap_combined_Lepr.png?v=4" width="33%" />

<img src="figures/umap_combined_Rspo3.png?v=4" width="33%" /><img src="figures/umap_combined_Sp7.png?v=4" width="33%" /><img src="figures/umap_combined_Tnfrsf11a.png?v=4" width="33%" />

<img src="figures/umap_combined_Kdr.png?v=4" width="33%" /><img src="figures/umap_combined_S100a8.png?v=4" width="33%" /><img src="figures/umap_combined_Abcc9.png?v=4" width="33%" />

<img src="figures/umap_combined_S100b.png?v=4" width="33%" /><img src="figures/umap_combined_S100a9.png?v=4" width="33%" /><img src="figures/umap_combined_Mbp.png?v=4" width="33%" />

<img src="figures/umap_combined_Vwf.png?v=4" width="33%" /><img src="figures/umap_combined_Myh11.png?v=4" width="33%" /><img src="figures/umap_combined_Krt5.png?v=4" width="33%" />

<img src="figures/umap_combined_Acp5.png?v=4" width="33%" /><img src="figures/umap_combined_Cdh1.png?v=4" width="33%" /><img src="figures/umap_combined_Mrc1.png?v=4" width="33%" />

<img src="figures/umap_combined_Frzb.png?v=4" width="33%" /><img src="figures/umap_combined_Kit.png?v=4" width="33%" /> <img src="figures/umap_combined_Adgre1.png?v=4" width="33%" />

<img src="figures/umap_combined_Comp.png?v=4" width="33%" /><img src="figures/umap_combined_Acan.png?v=4" width="33%" /> <img src="figures/umap_combined_Krt17.png?v=4" width="33%" />

<img src="figures/umap_combined_Hpgds.png?v=4" width="33%" /><img src="figures/umap_combined_Kcnj8.png?v=4" width="33%" /><img src="figures/umap_combined_Dmp1.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd28.png?v=4" width="33%" /><img src="figures/umap_combined_Csf3r.png?v=4" width="33%" /><img src="figures/umap_combined_Flt4.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd4.png?v=4" width="33%" /><img src="figures/umap_combined_Hdc.png?v=4" width="33%" /><img src="figures/umap_combined_Lgr6.png?v=4" width="33%" />

<img src="figures/umap_combined_Epcam.png?v=4" width="33%" /><img src="figures/umap_combined_Pi16.png?v=4" width="33%" /><img src="figures/umap_combined_Cilp2.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd163.png?v=4" width="33%" /><img src="figures/umap_combined_Dsg3.png?v=4" width="33%" /><img src="figures/umap_combined_Col9a1.png?v=4" width="33%" />

<img src="figures/umap_combined_Plp1.png?v=4" width="33%" /><img src="figures/umap_combined_Des.png?v=4" width="33%" /><img src="figures/umap_combined_Dcstamp.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd3e.png?v=4" width="33%" /><img src="figures/umap_combined_Col2a1.png?v=4" width="33%" /><img src="figures/umap_combined_Oscar.png?v=4" width="33%" />

<img src="figures/umap_combined_Lyve1.png?v=4" width="33%" /><img src="figures/umap_combined_Ighm.png?v=4" width="33%" /><img src="figures/umap_combined_Cd3d.png?v=4" width="33%" />

<img src="figures/umap_combined_Cnn1.png?v=4" width="33%" /><img src="figures/umap_combined_Aqp5.png?v=4" width="33%" /><img src="figures/umap_combined_Scnn1a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ngfr.png?v=4" width="33%" /><img src="figures/umap_combined_Prox1.png?v=4" width="33%" /><img src="figures/umap_combined_Calcr.png?v=4" width="33%" />

<img src="figures/umap_combined_Sp8.png?v=4" width="33%" /><img src="figures/umap_combined_Sp6.png?v=4" width="33%" /><img src="figures/umap_combined_Krt19.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1g.png?v=4" width="33%" /><img src="figures/umap_combined_Ms4a1.png?v=4" width="33%" /><img src="figures/umap_combined_Pax5.png?v=4" width="33%" />

<img src="figures/umap_combined_Sox10.png?v=4" width="33%" /><img src="figures/umap_combined_Mpo.png?v=4" width="33%" /><img src="figures/umap_combined_Cd22.png?v=4" width="33%" />

<img src="figures/umap_combined_Ucma.png?v=4" width="33%" /><img src="figures/umap_combined_Cd8a.png?v=4" width="33%" /><img src="figures/umap_combined_Matn1.png?v=4" width="33%" />

<img src="figures/umap_combined_Tpsb2.png?v=4" width="33%" /><img src="figures/umap_combined_Krt7.png?v=4" width="33%" /><img src="figures/umap_combined_Ccr3.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1b.png?v=4" width="33%" /><img src="figures/umap_combined_Ms4a2.png?v=4" width="33%" /><img src="figures/umap_combined_Krt18.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt8.png?v=4" width="33%" /><img src="figures/umap_combined_Cd19.png?v=4" width="33%" /><img src="figures/umap_combined_Foxa1.png?v=4" width="33%" />

<img src="figures/umap_combined_Ly6g.png?v=4" width="33%" /><img src="figures/umap_combined_Tpsab1.png?v=4" width="33%" /><img src="figures/umap_combined_Gdf5.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd200r3.png?v=4" width="33%" /><img src="figures/umap_combined_Mcpt8.png?v=4" width="33%" /><img src="figures/umap_combined_Cpa3.png?v=4" width="33%" />

<img src="figures/umap_combined_Bglap.png?v=4" width="33%" />

### Double checking cluster 5 

From the dotplot and feature plots,  cluster 5 does not show strong evidence of two distinct co-expressed gene programs.

![](figures/umap_doublet_score.png?v=1) 

From the dotplot, cluster 5 doesn't seem to be higlhy expressed in 2 distinct programs 
## doublet_score Summary Statistics

| Statistic | Value |
|----------|-------|
| Count    | 41304 |
| Mean     | 0.041292 |
| Std      | 0.045593 |
| Min      | 0.000852 |
| 25%      | 0.018277 |
| 50%      | 0.030149 |
| 75%      | 0.047481 |
| Max      | 0.484083 |


### Preliminary Annotations 

![](figures/umap_combined_celltypeON.png?v=8) 

![](figures/umap_combined_celltype.png?v=8)

### Now lets subset Fibroblast from GSE135985 and Osteosarcoma from AXT 

## Celltype Composition by Sample

| celltype      | NonReg_14DPA | Reg  | Reg_14DPA | Uninjured1 | Uninjured2 | nonReg |
|----------------|--------------|------|------------|-------------|-------------|---------|
| Fibroblast     | 634          | 0    | 433        | 1564        | 2521        | 0       |
| Osteosarcoma   | 0            | 1116 | 0          | 0           | 0           | 9075    |



