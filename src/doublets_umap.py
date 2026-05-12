import scanpy as sc
import matplotlib.pyplot as plt
import argparse

# -----------------------
# ARGPARSE INPUT
# -----------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
args = parser.parse_args()

# -----------------------
# LOAD DATA
# -----------------------
adata = sc.read_h5ad(args.input)

# -----------------------
# UMAP PLOT (continuous doublet score)
# -----------------------
sc.pl.umap(
    adata,
    color="doublet_score",   # continuous variable
    cmap="viridis",          # smooth gradient
    vmin=0,
    vmax=1,                  # Scrublet scores usually 0–1
    size=20,
    show=False
)

# -----------------------
# SAVE FIGURE
# -----------------------
plt.savefig(
   "figures/umap_doublet_score.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print((adata.obs["doublet_score"] > 0.5).sum())

print(adata.obs["doublet_score"].head())
print(adata.obs["doublet_score"].describe())

