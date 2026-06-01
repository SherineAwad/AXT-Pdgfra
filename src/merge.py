import scanpy as sc
import argparse
import os

# -------------------------
# Args
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input1', required=True)
parser.add_argument('--input2', required=True)
parser.add_argument('--output', required=True)

args = parser.parse_args()

# -------------------------
# Load datasets
# -------------------------
adata1 = sc.read_h5ad(args.input1)
adata2 = sc.read_h5ad(args.input2)

print(f"Input1: {adata1.n_obs} cells × {adata1.n_vars} genes")
print(f"Input2: {adata2.n_obs} cells × {adata2.n_vars} genes")

# -------------------------
# Ensure gene names are unique
# -------------------------
adata1.var_names_make_unique()
adata2.var_names_make_unique()

# -------------------------
# Add batch labels BEFORE merging
# -------------------------
batch1 = os.path.splitext(os.path.basename(args.input1))[0]
batch2 = os.path.splitext(os.path.basename(args.input2))[0]


batch1 = batch1.replace("_doublets", "")
batch2 = batch2.replace("_doublets", "")

adata1.obs["batch"] = batch1
adata2.obs["batch"] = batch2

# -------------------------
# Merge with OUTER join (keeps all genes)
# -------------------------
adata = sc.concat(
    [adata1, adata2],
    join="outer",
    fill_value=0,
    label="batch",
    keys=[batch1, batch2],
    index_unique="-"
)

print(f"Merged: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  (This includes all genes from both datasets)")

# -------------------------
# Optional: Verify no genes were lost
# -------------------------
original_total = len(set(adata1.var_names) | set(adata2.var_names))
print(f"  Expected unique genes: {original_total}")
print(f"  Actual genes in merge: {adata.n_vars}")

# -------------------------
# Save
# -------------------------
adata.write(args.output, compression="gzip")

print(f"SAVED → {args.output}")
