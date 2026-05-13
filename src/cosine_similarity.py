#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse

warnings.filterwarnings('ignore')


def compute_state_pseudobulk(adata, celltype_col, sample_col):
    """
    Each state = (celltype, sample)
    """
    states = adata.obs[celltype_col].astype(str) + "_" + adata.obs[sample_col].astype(str)
    adata.obs["state"] = states

    profiles = {}

    for state in adata.obs["state"].unique():
        mask = adata.obs["state"] == state
        X = adata[mask].X

        if sparse.issparse(X):
            profiles[state] = np.array(X.mean(axis=0)).flatten()
        else:
            profiles[state] = X.mean(axis=0)

    return profiles


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)

    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    celltype_col = "celltype"
    sample_col = "sample"

    # -------------------------------------------------------
    # HVG selection on full dataset
    # -------------------------------------------------------
    print("Selecting HVGs...")

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat"
    )

    genes = adata.var_names[adata.var["highly_variable"]]
    adata = adata[:, genes].copy()

    # -------------------------------------------------------
    # PSEUDOBULK per (celltype, sample)
    # -------------------------------------------------------
    print("Computing state-level pseudobulk...")

    profiles = compute_state_pseudobulk(adata, celltype_col, sample_col)

    states = list(profiles.keys())
    mat = np.vstack([profiles[s] for s in states])

    # -------------------------------------------------------
    # COSINE SIMILARITY (FULL MATRIX)
    # -------------------------------------------------------
    print("Computing similarity matrix...")

    sim_matrix = cosine_similarity(mat, mat)

    sim_df = pd.DataFrame(
        sim_matrix,
        index=states,
        columns=states
    )

    # -------------------------------------------------------
    # SAVE
    # -------------------------------------------------------
    csv_path = f"{args.prefix}_cosine_similarity.csv"
    sim_df.to_csv(csv_path)

    print(f"Saved CSV: {csv_path}")

    # -------------------------------------------------------
    # PLOT
    # -------------------------------------------------------
    print("Plotting...")

    plt.figure(figsize=(
        max(10, len(states) * 0.6),
        max(10, len(states) * 0.6)
    ))

    plt.imshow(sim_df.values, cmap="viridis", aspect="auto")
    plt.colorbar(label="Cosine similarity")
    data = sim_df.values

    for i in range(data.shape[0]):
       for j in range(data.shape[1]):
          plt.text(
              j, i,
              f"{data[i, j]:.2f}",
              ha='center',
              va='center',
              color='black',
              fontsize=5  # adjust if too dense
           )
    plt.xticks(range(len(states)), states, rotation=90, fontsize=6)
    plt.yticks(range(len(states)), states, fontsize=6)

    plt.title("Celltype × Sample State Similarity Map")

    plt.tight_layout()

    os.makedirs("figures", exist_ok=True)

    plt.savefig(
        f"figures/{args.prefix}_cosine_similarity.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Done.")


if __name__ == "__main__":
    main()
