#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.linear_model import LinearRegression

# -------------------------
# ARGPARSE
# -------------------------
parser = argparse.ArgumentParser(description="Quantify contribution of EGFP to PCA")
parser.add_argument("--input", required=True, help="Input .h5ad file")
parser.add_argument("--prefix", required=True, help="Prefix for output CSV and figures")
args = parser.parse_args()

# -------------------------
# LOAD DATA
# -------------------------
adata = sc.read_h5ad(args.input)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

outdir = os.path.abspath("figures")
os.makedirs(outdir, exist_ok=True)

# -------------------------
# PCA
# -------------------------
if "X_pca" not in adata.obsm:
    sc.pp.pca(adata, n_comps=20)
pc1 = adata.obsm["X_pca"][:, 0]
pc2 = adata.obsm["X_pca"][:, 1]

# -------------------------
# EGFP EXPRESSION
# -------------------------
if "EGFP" not in adata.var_names:
    raise ValueError("EGFP not found in adata.var_names")

egfp = adata[:, "EGFP"].X
if hasattr(egfp, "toarray"):
    egfp = egfp.toarray().ravel()
else:
    egfp = np.ravel(egfp)

# -------------------------
# LINEAR REGRESSION (VARIANCE EXPLAINED)
# -------------------------
def reg_r2(y, X):
    """Return R² of linear regression of y on X"""
    if np.std(X) == 0:
        return np.nan
    model = LinearRegression().fit(X.reshape(-1, 1), y)
    return model.score(X.reshape(-1, 1), y)

r2_pc1 = reg_r2(pc1, egfp)
r2_pc2 = reg_r2(pc2, egfp)

print(f"EGFP -> PC1 R²: {r2_pc1:.3f}")
print(f"EGFP -> PC2 R²: {r2_pc2:.3f}")

# -------------------------
# SAVE CSV
# -------------------------
df = pd.DataFrame({
    "PC": ["PC1", "PC2"],
    "R2_EGFP": [r2_pc1, r2_pc2]
})
csv_path = os.path.join(outdir, f"{args.prefix}_egfp_pca_r2.csv")
df.to_csv(csv_path, index=False)
print("Saved CSV:", csv_path)

# -------------------------
# OPTIONAL FIGURE: Scatter + regression line
# -------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for i, pc in enumerate([pc1, pc2]):
    ax[i].scatter(egfp, pc, s=5, alpha=0.6)
    # regression line
    if np.std(egfp) > 0:
        coef = np.polyfit(egfp, pc, 1)
        xs = np.array([egfp.min(), egfp.max()])
        ys = coef[0] * xs + coef[1]
        ax[i].plot(xs, ys, color='red', lw=2)
    ax[i].set_xlabel("EGFP expression")
    ax[i].set_ylabel(f"PC{i+1} score")
    ax[i].set_title(f"PC{i+1} vs EGFP (R²={df.loc[i,'R2_EGFP']:.3f})")
plt.tight_layout()
fig_path = os.path.join(outdir, f"{args.prefix}_egfp_pc_r2.png")
plt.savefig(fig_path, dpi=120, bbox_inches="tight")
plt.close()
print("Saved figure:", fig_path)
