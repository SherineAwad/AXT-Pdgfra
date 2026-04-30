import scanpy as sc

# Read each sample
adata1 = sc.read_10x_mtx("8_Processed_Uninjured2", var_names='gene_symbols')
adata2 = sc.read_10x_mtx("7_Processed_Uninjured1", var_names='gene_symbols')
adata3 = sc.read_10x_mtx("3_Processed_14DPA", var_names='gene_symbols')
adata4 = sc.read_10x_mtx("11_Processed_NonRegen_14DPA", var_names='gene_symbols')

# Add sample labels (IMPORTANT to track origin)
adata1.obs['sample'] = 'Uninjured2'
adata2.obs['sample'] = 'Uninjured1'
adata3.obs['sample'] = 'Reg_14DPA'
adata4.obs['sample'] = 'NonReg_14DPA'

# Make gene names unique (avoids merge issues)
for ad in [adata1, adata2, adata3, adata4]:
    ad.var_names_make_unique()

# Merge all datasets
adata = adata1.concatenate(
    adata2, adata3, adata4,
    batch_key='batch',
    batch_categories=['Uninjured2', 'Uninjured1', 'Reg_14DPA', 'NonReg_14DPA']
)

# Save to h5ad
adata.write("GSE135985.h5ad")

print("Done")
