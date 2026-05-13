#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import warnings

warnings.simplefilter("ignore", RuntimeWarning)


# ----------------------------
# GET MATRIX
# ----------------------------
def get_X(adata):
    X = adata.X
    return X.toarray() if issparse(X) else X


# ----------------------------
# MAIN
# ----------------------------
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--max_cells", type=int, default=10000)

    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    assert "celltype" in adata.obs
    assert "sample" in adata.obs

    # ----------------------------
    # BUILD STATES (celltype × sample)
    # ----------------------------
    adata.obs["state"] = (
        adata.obs["celltype"].astype(str)
        + "_"
        + adata.obs["sample"].astype(str)
    )

    states = sorted(adata.obs["state"].unique())
    print(f"Found {len(states)} states")

    rng = np.random.default_rng(42)

    # ----------------------------
    # COLLECT CELLS PER STATE
    # ----------------------------
    state_data = {}

    all_X = []

    for st in states:
        sub = adata[adata.obs["state"] == st]

        if sub.n_obs == 0:
            continue

        if sub.n_obs > args.max_cells:
            idx = rng.choice(sub.n_obs, args.max_cells, replace=False)
            sub = sub[idx]

        X = get_X(sub)
        state_data[st] = X
        all_X.append(X)

    X_all = np.vstack(all_X)

    print(f"Total cells used: {X_all.shape}")

    # ----------------------------
    # PCA FIT
    # ----------------------------
    print("Fitting PCA...")

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_all)

    pca = PCA(n_components=min(args.n_pcs, X_scaled.shape[1], X_scaled.shape[0]))
    pca.fit(X_scaled)

    print(f"PCA components: {pca.n_components_}")

    # ----------------------------
    # PROJECT STATES
    # ----------------------------
    proj = {}

    start = 0

    for st in states:
        X = state_data[st]
        n = X.shape[0]

        X_scaled = scaler.transform(X)
        proj[st] = pca.transform(X_scaled)

        start += n

    # ----------------------------
    # WASSERSTEIN DISTANCE MATRIX
    # ----------------------------
    n_st = len(states)
    dist = np.zeros((n_st, n_st))

    print("Computing Wasserstein distances...")

    for i, s1 in enumerate(states):
        for j, s2 in enumerate(states):

            A = proj[s1]
            B = proj[s2]

            n_pcs = min(A.shape[1], B.shape[1])

            dists = [
                wasserstein_distance(A[:, k], B[:, k])
                for k in range(n_pcs)
            ]

            dist[i, j] = np.mean(dists)

        print(f"Done: {s1}")

    # ----------------------------
    # SYMMETRIZE
    # ----------------------------
    dist = (dist + dist.T) / 2

    # ----------------------------
    # CONVERT TO SIMILARITY
    # ----------------------------
    max_d = dist.max()
    sim = 1 - dist / max_d if max_d > 0 else np.ones_like(dist)

    # ----------------------------
    # PLOT
    # ----------------------------
    os.makedirs("figures", exist_ok=True)

    plt.figure(figsize=(max(10, n_st * 0.6), max(10, n_st * 0.6)))

    plt.imshow(sim, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    plt.colorbar(label="PCA + Wasserstein similarity")

    plt.xticks(range(n_st), states, rotation=90, fontsize=6)
    plt.yticks(range(n_st), states, fontsize=6)

    # annotations
    for i in range(n_st):
        for j in range(n_st):
            plt.text(
                j, i,
                f"{sim[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=5,
                color="black"
            )

    plt.title(f"PCA + Wasserstein State Similarity\n{args.prefix}")

    plt.tight_layout()

    out = f"figures/{args.prefix}_pca_wasserstein_similarity.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
