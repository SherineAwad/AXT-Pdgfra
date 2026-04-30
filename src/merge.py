import scanpy as sc
import argparse

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
# Ensure gene alignment
# -------------------------
adata1.var_names_make_unique()
adata2.var_names_make_unique()

# keep only shared genes (safe merge)
shared_genes = adata1.var_names.intersection(adata2.var_names)

adata1 = adata1[:, shared_genes].copy()
adata2 = adata2[:, shared_genes].copy()

print(f"Shared genes: {len(shared_genes)}")

# -------------------------
# Add batch labels
# -------------------------
adata1.obs["batch"] = "input1"
adata2.obs["batch"] = "input2"

# -------------------------
# Merge
# -------------------------
adata = sc.concat(
    [adata1, adata2],
    join="inner",
    label="batch",
    keys=["input1", "input2"],
    index_unique="-"
)

print(f"Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# -------------------------
# Save
# -------------------------
adata.write(args.output, compression="gzip")

print(f"SAVED → {args.output}")
