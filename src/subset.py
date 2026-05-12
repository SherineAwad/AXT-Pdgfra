import scanpy as sc
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True, help="Input h5ad file")
parser.add_argument("--output", required=True, help="Output h5ad file")

parser.add_argument(
    "--subset",
    required=True,
    help=(
        "Subset definition in format: "
        "'T cells:A,B;B cells:C,D,E,F'"
    )
)

parser.add_argument(
    "--prefix",
    default="output",
    help="Prefix for saved figures"
)

args = parser.parse_args()

# -------------------------
# Read data
# -------------------------
adata = sc.read_h5ad(args.input)

# -------------------------
# Validation
# -------------------------
if "celltype" not in adata.obs.columns:
    raise ValueError("celltype column not found in adata.obs")

if "sample" not in adata.obs.columns:
    raise ValueError("sample column not found in adata.obs")

# -------------------------
# Parse subset rules
# -------------------------
subset_rules = {}

try:
    rules = args.subset.split(";")

    for rule in rules:
        celltype, samples = rule.split(":")
        celltype = celltype.strip()

        sample_list = [s.strip() for s in samples.split(",")]

        subset_rules[celltype] = sample_list

except Exception:
    raise ValueError(
        "Invalid --subset format. Example:\n"
        "'T cells:A,B;B cells:C,D,E,F'"
    )

# -------------------------
# Validate celltypes/samples
# -------------------------
available_celltypes = adata.obs["celltype"].unique()
available_samples = adata.obs["sample"].unique()

for celltype, samples in subset_rules.items():

    # validate celltype
    if celltype not in available_celltypes:
        raise ValueError(
            f"Celltype '{celltype}' not found.\n"
            f"Available: {available_celltypes}"
        )

    # validate samples
    invalid_samples = [
        s for s in samples
        if s not in available_samples
    ]

    if invalid_samples:
        raise ValueError(
            f"Sample(s) {invalid_samples} not found.\n"
            f"Available: {available_samples}"
        )

    # validate celltype exists in requested samples
    existing_samples = adata.obs.loc[
        adata.obs["celltype"] == celltype,
        "sample"
    ].unique()

    missing = [
        s for s in samples
        if s not in existing_samples
    ]

    if missing:
        raise ValueError(
            f"Celltype '{celltype}' not found in sample(s): {missing}"
        )

# -------------------------
# Build mask
# -------------------------
mask = np.zeros(adata.n_obs, dtype=bool)

for celltype, samples in subset_rules.items():

    current_mask = (
        (adata.obs["celltype"] == celltype)
        &
        (adata.obs["sample"].isin(samples))
    )

    mask |= current_mask

# -------------------------
# Subset data
# -------------------------
adata_subset = adata[mask].copy()

print(
    f"Subsetted data: "
    f"{adata_subset.n_obs} cells, "
    f"{adata_subset.n_vars} genes"
)

# -------------------------
# Create output dir
# -------------------------
os.makedirs("figures", exist_ok=True)

# -------------------------
# GLOBAL UMAP
# -------------------------
fig, ax = plt.subplots(figsize=(6, 5))

sc.pl.umap(
    adata_subset,
    color="sample",
    size=20,
    ax=ax,
    show=False
)

plt.tight_layout()
plt.savefig(
    f"figures/{args.prefix}_umap.png",
    dpi=300
)
plt.close()

# -------------------------
# PER-SAMPLE UMAP
# -------------------------
samples = adata_subset.obs["sample"].unique()

fig, axes = plt.subplots(
    1,
    len(samples),
    figsize=(5 * len(samples), 5)
)

# handle single sample case
if len(samples) == 1:
    axes = [axes]

for i, s in enumerate(samples):

    sc.pl.umap(
        adata_subset[
            adata_subset.obs["sample"] == s
        ],
        color="sample",
        size=20,
        title=f"{s}",
        ax=axes[i],
        show=False
    )

plt.tight_layout()

plt.savefig(
    f"figures/{args.prefix}_perSample_umap.png",
    dpi=300
)

plt.close()

# -------------------------
# Save subset
# -------------------------
adata_subset.write_h5ad(args.output)

print(f"Saved: {args.output}")

# -------------------------
# Summary
# -------------------------
print("\nSummary:")
print(f"  Cells: {adata_subset.n_obs}")
print(f"  Genes: {adata_subset.n_vars}")

print("\nCelltype x Sample breakdown:")

summary_table = (
    adata_subset.obs
    .groupby(["celltype", "sample"])
    .size()
    .unstack(fill_value=0)
)

print(summary_table)


