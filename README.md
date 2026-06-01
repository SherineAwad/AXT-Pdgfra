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
![](figures/violin_GSE135985_preQC.png?v=2)

#### Post filtering 
![](figures/violin_GSE135985_AfterQC.png?v=2)

### Doublet detection 

![](figures/GSE135985_scrublet_scores.png?v=2)

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


![](figures/umap_combined_umap.png?v=5)

<img src="figures/umap_combined_nonReg.png?v=5" width="33%" /><img src="figures/umap_combined_Reg_14DPA.png?v=4" width="33%" /><img src="figures/umap_combined_Uninjured1.png?v=4" width="33%" />

<img src="figures/umap_combined_Uninjured2.png?v=5" width="33%" /><img src="figures/umap_combined_Reg.png?v=4" width="33%" /><img src="figures/umap_combined_NonReg_14DPA.png?v=4" width="33%" />

## Use run_harmony for batch effect 

##  🟢 🟢 🟢 Harmony using All Genes ~30k genes

![](figures/umap_harmony_noHVG_harmony_batch.png?v=2) 


### Clustering 

![](figures/umap_combined_leiden.png?v=5)

<img src="figures/violin_combined_QC_n_genes_by_counts.png?v=5" width="33%" /><img src="figures/violin_combined_QC_total_counts.png?v=5" width="33%" /><img src="figures/violin_combined_QC_pct_counts_mt.png?v=5" width="33%" />

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

<img src="figures/umap_combined_Nfatc1.png?v=5" width="33%" /><img src="figures/umap_combined_Col1a1.png?v=4" width="33%" /><img src="figures/umap_combined_Pecam1.png?v=4" width="33%" />

<img src="figures/umap_combined_Prrx1.png?v=5" width="33%" /><img src="figures/umap_combined_Spp1.png?v=4" width="33%" /><img src="figures/umap_combined_Cd44.png?v=4" width="33%" />

<img src="figures/umap_combined_Dcn.png?v=5" width="33%" /><img src="figures/umap_combined_Runx2.png?v=4" width="33%" /><img src="figures/umap_combined_Postn.png?v=4" width="33%" />

<img src="figures/umap_combined_Acta2.png?v=5" width="33%" /><img src="figures/umap_combined_Pdgfrb.png?v=4" width="33%" /><img src="figures/umap_combined_Ibsp.png?v=4" width="33%" />

<img src="figures/umap_combined_Pmp22.png?v=5" width="33%" /><img src="figures/umap_combined_Cspg4.png?v=4" width="33%" /><img src="figures/umap_combined_Cd79a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ctsk.png?v=5" width="33%" /><img src="figures/umap_combined_Rgs5.png?v=4" width="33%" /><img src="figures/umap_combined_Pdgfra.png?v=4" width="33%" />

<img src="figures/umap_combined_EGFP.png?v=5" width="33%" /><img src="figures/umap_combined_Cxcl12.png?v=4" width="33%" /><img src="figures/umap_combined_Mmp13.png?v=4" width="33%" />

<img src="figures/umap_combined_Col23a1.png?v=5" width="33%" /><img src="figures/umap_combined_Alpl.png?v=4" width="33%" /><img src="figures/umap_combined_Erg.png?v=4" width="33%" />

<img src="figures/umap_combined_Flt1.png?v=5" width="33%" /><img src="figures/umap_combined_Sox9.png?v=4" width="33%" /><img src="figures/umap_combined_Prg4.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt14.png?v=5" width="33%" /><img src="figures/umap_combined_Tagln.png?v=4" width="33%" /><img src="figures/umap_combined_Cd34.png?v=4" width="33%" />

<img src="figures/umap_combined_Esam.png?v=5" width="33%" /><img src="figures/umap_combined_Pdpn.png?v=4" width="33%" /><img src="figures/umap_combined_Cd68.png?v=4" width="33%" />

<img src="figures/umap_combined_Cdh5.png?v=5" width="33%" /><img src="figures/umap_combined_Emcn.png?v=4" width="33%" /><img src="figures/umap_combined_Plvap.png?v=4" width="33%" />

<img src="figures/umap_combined_Nes.png?v=5" width="33%" /><img src="figures/umap_combined_Csf1r.png?v=4" width="33%" /><img src="figures/umap_combined_Lepr.png?v=4" width="33%" />

<img src="figures/umap_combined_Rspo3.png?v=5" width="33%" /><img src="figures/umap_combined_Sp7.png?v=4" width="33%" /><img src="figures/umap_combined_Tnfrsf11a.png?v=4" width="33%" />

<img src="figures/umap_combined_Kdr.png?v=5" width="33%" /><img src="figures/umap_combined_S100a8.png?v=4" width="33%" /><img src="figures/umap_combined_Abcc9.png?v=4" width="33%" />

<img src="figures/umap_combined_S100b.png?v=5" width="33%" /><img src="figures/umap_combined_S100a9.png?v=4" width="33%" /><img src="figures/umap_combined_Mbp.png?v=4" width="33%" />

<img src="figures/umap_combined_Vwf.png?v=5" width="33%" /><img src="figures/umap_combined_Myh11.png?v=4" width="33%" /><img src="figures/umap_combined_Krt5.png?v=4" width="33%" />

<img src="figures/umap_combined_Acp5.png?v=5" width="33%" /><img src="figures/umap_combined_Cdh1.png?v=4" width="33%" /><img src="figures/umap_combined_Mrc1.png?v=4" width="33%" />

<img src="figures/umap_combined_Frzb.png?v=5" width="33%" /><img src="figures/umap_combined_Kit.png?v=4" width="33%" /> <img src="figures/umap_combined_Adgre1.png?v=4" width="33%" />

<img src="figures/umap_combined_Comp.png?v=5" width="33%" /><img src="figures/umap_combined_Acan.png?v=4" width="33%" /> <img src="figures/umap_combined_Krt17.png?v=4" width="33%" />

<img src="figures/umap_combined_Hpgds.png?v=5" width="33%" /><img src="figures/umap_combined_Kcnj8.png?v=4" width="33%" /><img src="figures/umap_combined_Dmp1.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd28.png?v=5" width="33%" /><img src="figures/umap_combined_Csf3r.png?v=4" width="33%" /><img src="figures/umap_combined_Flt4.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd4.png?v=5" width="33%" /><img src="figures/umap_combined_Hdc.png?v=4" width="33%" /><img src="figures/umap_combined_Lgr6.png?v=4" width="33%" />

<img src="figures/umap_combined_Epcam.png?v=5" width="33%" /><img src="figures/umap_combined_Pi16.png?v=4" width="33%" /><img src="figures/umap_combined_Cilp2.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd163.png?v=5" width="33%" /><img src="figures/umap_combined_Dsg3.png?v=4" width="33%" /><img src="figures/umap_combined_Col9a1.png?v=4" width="33%" />

<img src="figures/umap_combined_Plp1.png?v=5" width="33%" /><img src="figures/umap_combined_Des.png?v=4" width="33%" /><img src="figures/umap_combined_Dcstamp.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd3e.png?v=5" width="33%" /><img src="figures/umap_combined_Col2a1.png?v=4" width="33%" /><img src="figures/umap_combined_Oscar.png?v=4" width="33%" />

<img src="figures/umap_combined_Lyve1.png?v=5" width="33%" /><img src="figures/umap_combined_Ighm.png?v=4" width="33%" /><img src="figures/umap_combined_Cd3d.png?v=4" width="33%" />

<img src="figures/umap_combined_Cnn1.png?v=5" width="33%" /><img src="figures/umap_combined_Aqp5.png?v=4" width="33%" /><img src="figures/umap_combined_Scnn1a.png?v=4" width="33%" />

<img src="figures/umap_combined_Ngfr.png?v=5" width="33%" /><img src="figures/umap_combined_Prox1.png?v=4" width="33%" /><img src="figures/umap_combined_Calcr.png?v=4" width="33%" />

<img src="figures/umap_combined_Sp8.png?v=5" width="33%" /><img src="figures/umap_combined_Sp6.png?v=4" width="33%" /><img src="figures/umap_combined_Krt19.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1g.png?v=5" width="33%" /><img src="figures/umap_combined_Ms4a1.png?v=4" width="33%" /><img src="figures/umap_combined_Pax5.png?v=4" width="33%" />

<img src="figures/umap_combined_Sox10.png?v=5" width="33%" /><img src="figures/umap_combined_Mpo.png?v=4" width="33%" /><img src="figures/umap_combined_Cd22.png?v=4" width="33%" />

<img src="figures/umap_combined_Ucma.png?v=5" width="33%" /><img src="figures/umap_combined_Cd8a.png?v=4" width="33%" /><img src="figures/umap_combined_Matn1.png?v=4" width="33%" />

<img src="figures/umap_combined_Tpsb2.png?v=5" width="33%" /><img src="figures/umap_combined_Krt7.png?v=4" width="33%" /><img src="figures/umap_combined_Ccr3.png?v=4" width="33%" />

<img src="figures/umap_combined_Scnn1b.png?v=5" width="33%" /><img src="figures/umap_combined_Ms4a2.png?v=4" width="33%" /><img src="figures/umap_combined_Krt18.png?v=4" width="33%" />

<img src="figures/umap_combined_Krt8.png?v=5" width="33%" /><img src="figures/umap_combined_Cd19.png?v=4" width="33%" /><img src="figures/umap_combined_Foxa1.png?v=4" width="33%" />

<img src="figures/umap_combined_Ly6g.png?v=5" width="33%" /><img src="figures/umap_combined_Tpsab1.png?v=4" width="33%" /><img src="figures/umap_combined_Gdf5.png?v=4" width="33%" />

<img src="figures/umap_combined_Cd200r3.png?v=5" width="33%" /><img src="figures/umap_combined_Mcpt8.png?v=4" width="33%" /><img src="figures/umap_combined_Cpa3.png?v=4" width="33%" />

<img src="figures/umap_combined_Bglap.png?v=5" width="33%" />

### Double checking cluster doublets 

From the dotplot and feature plots, 

![](figures/umap_doublet_score.png?v=2) 



### Preliminary Annotations 

![](figures/umap_combined_celltypeON.png?v=11) 

![](figures/umap_combined_celltype.png?v=11)


# From here requires re-running -- ocming soon 
### Now lets subset Fibroblast from GSE135985 and Osteosarcoma from AXT 

## Celltype Composition by Sample

| celltype      | Non_Regen_14_DPA | Reg | Regen_14_DPA | Uninjured1 | Uninjured2 | nonReg |
|---------------|------------------|-----|--------------|-------------|-------------|--------|
| Fibroblast    | 744              | 0   | 638          | 1663        | 3011        | 0      |
| Osteosarcoma  | 0                | 1108| 0            | 0           | 0           | 8940   |


![](figures/Fibroblast_Osteosarcoma_umap.png?v=2)

![](figures/Fibroblast_Osteosarcoma_perSample_umap.png?v=2)

## Celltype similarity 

### Method 1: Pearson Correlation 

> whether genes vary together above or below their average in both states

- centers each vector (subtracts mean)
- captures linear co-variation
- sensitive to relative up/down shifts

Interpretation:
> similarity of co-activation patterns relative to mean expression
> “Do genes go up/down together relative to average?”

![](figures/Fibroblast_Osteosarcoma_pearson_hvg2000_similarity.png?v=1)

### Method 2: Spearman Correlation 

> whether two states agree on the ranking of genes

- converts values into ranks first
- ignores magnitude completely
- compares ordering of genes

Interpretation:
> similarity of gene importance or ranking
> “Do they agree on which genes are most important?”


![](figures/Fibroblast_Osteosarcoma_spearman_hvg2000_similarity.png?v=1)

### Method 3: Cosine similarity 

##### 1. What is a “vector” in this analysis?

Each biological state is defined as:

(celltype × sample)

Example:
- Osteosarcoma_Reg
- Fibroblast_Injured1

Each state is converted into a gene expression vector.

If we have 10 genes:

A = [
  mean expression of gene1 in Osteosarcoma_Reg,
  mean expression of gene2 in Osteosarcoma_Reg,
  mean expression of gene3 in Osteosarcoma_Reg,
  ...
  mean expression of gene10 in Osteosarcoma_Reg
]

This is called a pseudobulk expression profile.

So:
- Vector A = one state (e.g. Osteosarcoma_Reg)
- Vector B = another state (e.g. Fibroblast_Injured1)

---

##### 2. What cosine similarity compares

Cosine similarity compares:

the pattern of gene expression between two states

It does NOT depend on:
- number of cells
- absolute expression scale
- sequencing depth

It only cares about:
whether genes go up and down together in both states

---

##### 3. Mathematical definition

Given two vectors A and B:

cos(A, B) = (A · B) / (||A|| × ||B||)

Where:

Dot product:
A · B = Σ (Ai × Bi)

Norm (length):
||A|| = sqrt(Σ Ai²)
||B|| = sqrt(Σ Bi²)

---

##### 4. What this means intuitively

Cosine similarity asks:

“Do these two gene expression profiles point in the same direction?”

- Value = 1 → identical pattern
- Value = 0 → unrelated patterns
- Value = -1 → opposite patterns (rare in expression data)

---

##### 5. Biological interpretation

Each vector is a cell state signature.

Cosine similarity tells us:

how similar two biological states are in their transcriptional programs; the alignment of gene expression programs between biological states, ignoring scale and focusing only on expression pattern shape

![](figures/Fibroblast_Osteosarcoma_cosine_hvg2000_similarity.png?v=1)

### Method 4: PCA + Wasserstein Distance 

#####  What is being compared?

Each data point is a **state**:

(celltype × sample)

Example states:
- Osteosarcoma_Reg
- Osteosarcoma_nonReg
- Fibroblast_14DPA
- Fibroblast_Uninjured1

Each state contains many single cells.

So instead of comparing single vectors, we compare:
> **distributions of cells**

---

##### Step 1: Represent cells in a shared space (PCA)

All cells from all states are combined and projected into PCA space.

Each cell becomes:

(x₁, x₂, x₃, ..., xₙ)

Where:
- each axis = principal component
- PCA captures main sources of variation in gene expression

So now:
> each state = a cloud of points in PCA space

---

##### Step 2: Each state becomes a distribution

For each state:

- collect all its cells
- represent them in PCA space

So instead of a single vector, you now have:

> a **distribution of points per state**

Example:
- Osteosarcoma_Reg → cloud A
- Fibroblast_Injured1 → cloud B

---

##### Step 3: Compare distributions (Wasserstein distance)

Wasserstein distance measures:

> how much “work” is needed to move one distribution into another

Intuition:
- if two clouds overlap → low distance
- if they are far apart → high distance

It compares:
> shape + spread + location of cell populations

---

##### Step 4: Per-PC comparison

Instead of comparing full high-dimensional distributions at once, your implementation:

- compares each PCA dimension separately
- computes Wasserstein distance per PC
- averages across PCs


---

##### Step 5: Convert distance → similarity

Since Wasserstein is a **distance (bigger = more different)**:

similarity = 1 - normalized_distance

So:
- 1 → identical distributions
- 0 → completely different distributions


![](figures/Fibroblast_Osteosarcom_pca_wasserstein_hvg2000.png?v=1)


## Method 5: PCA-based Maximum Mean Discrepancy (MMD)

##### PCA step (same as Wasserstein)
- one global PCA is fitted using all cells from all celltypes and samples
- all cells are projected into this shared PCA space
- each state = (celltype × sample) becomes a cloud of points in PCA space

---

##### MMD step

Instead of measuring transport distance between two clouds (as in Wasserstein), MMD:

- takes two states A and B (e.g., Osteosarcoma_Reg vs Fibroblast_Uninjured1)
- compares the statistical difference between the two clouds in PCA space using a kernel function
- produces a single value that reflects how different the two distributions are

![](figures/Fibroblast_Osteosarcom_pca_mmd_hvg2000.png?v=1)

## Method 6: Python Optimal Transport 

POT answers one question: **What is the minimum cost to transform one distribution of cells into another distribution of cells?**

##### How It Works

**Step 1: PCA**
- Reduce gene expression data to ~50 principal components
- Each cell becomes a point in PCA space

**Step 2: Build Cost Matrix**
- Calculate distance between every cell in group A and every cell in group B
- This matrix (M) represents how "expensive" it is to move mass from any A cell to any B cell

**Step 3: Sinkhorn Algorithm**
- Finds the optimal transport plan (T)
- T tells you how much mass from each A cell should go to each B cell
- Minimizes total transport cost

**Step 4: Calculate Similarity**
- similarity = 1 - (total transport cost)
- High similarity (close to 1) = distributions overlap well
- Low similarity (close to 0) = distributions are very different

##### Why POT for Cell Types

- Different number of cells per sample? OT handles unequal group sizes naturally
- Cell states exist on a continuum? Compares entire distributions, not just averages
- Two populations can have same mean but different structure? Captures differences in spread, shape, and density

##### Analogy

Imagine two clouds of points in PCA space. If the clouds heavily overlap, low transport cost and cells are similar. If the clouds are far apart, high transport cost and cells are different. POT measures the "work" needed to morph one cloud into the other.

![](figures/Fibroblast_Osteosarcoma_pot_hvg2000_matrix.png?v=1)

| Aspect | Cosine | Pearson | Spearman | PCA + MMD | PCA + Wasserstein | PCA + POT (Sinkhorn) |
|--------|--------|---------|----------|-----------|-------------------|----------------------|
| What it reflects | Similarity of which genes are high vs low | Similarity of how much genes go up and down together | Similarity of which genes are highest and lowest ranked | Whether two cell populations have the same cellular composition | Whether two cell populations have the same cellular composition | Whether two cell populations have the same cellular composition |
| What it tells you | Do the same genes have high expression in both groups? | When one gene goes up, does the same gene go up in the other group? | Do both groups rank genes from most expressed to least expressed in the same order? | Do both groups contain the same cell types and states in the same proportions? | Do both groups contain the same cell types and states in the same proportions? | Do both groups contain the same cell types and states in the same proportions? |


### 🚨🚨🚨 Method 7 SCOT+: SCOT+ is missing -- software issue -- contacting author 





