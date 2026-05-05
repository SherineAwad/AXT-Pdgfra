#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    # -------------------------
    # LOAD
    # -------------------------
    adata = sc.read(args.input)

    # 🔥 KEEP FULL COPY BEFORE ANYTHING
    adata_full = adata.copy()

    outdir = os.path.abspath("figures")
    os.makedirs(outdir, exist_ok=True)
    print("Saving to:", outdir)

    # -------------------------
    # SAFE column handling
    # -------------------------
    var_df = adata.var.copy()

    # Try different possible column names
    if "means" in var_df.columns:
        mean = var_df["means"]
    elif "mean" in var_df.columns:
        mean = var_df["mean"]
    else:
        # Compute mean if not present
        mean = np.array(adata.X.mean(axis=0)).flatten()
        var_df["means"] = mean

    if "variances_norm" in var_df.columns:
        var = var_df["variances_norm"]
    elif "variances" in var_df.columns:
        var = var_df["variances"]
    elif "dispersions_norm" in var_df.columns:
        var = var_df["dispersions_norm"]
    else:
        # Compute variance for sparse matrix
        X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        var = np.var(X_dense, axis=0)
        var_df["variances"] = var

    if "highly_variable" in var_df.columns:
        hvg_flag = var_df["highly_variable"]
    else:
        # Mark all genes as non-HVG if column missing
        hvg_flag = pd.Series([False] * len(var_df), index=var_df.index)

    # -------------------------
    # Save HVG list
    # -------------------------
    hvg = var_df[hvg_flag].copy()
    hvg = hvg.reset_index().rename(columns={"index": "gene"})

    csv_path = os.path.join(outdir, f"{args.prefix}_hvg.csv")
    hvg.to_csv(csv_path, index=False)

    # -------------------------
    # Build dataframe
    # -------------------------
    df = pd.DataFrame({
        "gene": adata.var_names,
        "mean": mean,
        "variance": var,
        "hvg": hvg_flag
    })

    # -------------------------
    # HVG plot
    # -------------------------
    plt.figure(figsize=(6, 5))
    plt.scatter(df["mean"], df["variance"], s=2, alpha=0.2)
    if df["hvg"].any():
        plt.scatter(df[df["hvg"]]["mean"], df[df["hvg"]]["variance"], s=5)
    plt.xlabel("Mean")
    plt.ylabel("Variance")
    plt.title("HVG selection")
    plt.savefig(os.path.join(outdir, f"{args.prefix}_hvg.png"), dpi=150)
    plt.close()

    # =====================================================
    # 🔥 EGFP FROM FULL DATA (THIS IS THE FIX)
    # =====================================================
    egfp_gene = None

    for g in adata_full.var_names:
        if "egfp" in g.lower() or g.lower() == "gfp":
            egfp_gene = g
            break

    if egfp_gene is None:
        print("❌ EGFP NOT FOUND IN FULL DATA")
        print("Example genes:", list(adata_full.var_names[:20]))
        return

    print(f"✅ Found EGFP: {egfp_gene}")

    egfp = adata_full[:, egfp_gene].X

    if hasattr(egfp, "toarray"):
        egfp = egfp.toarray().flatten()
    else:
        egfp = np.array(egfp).flatten()

    print("EGFP mean:", np.mean(egfp))
    print("EGFP std:", np.std(egfp))
    print("Non-zero cells:", np.sum(egfp > 0), "/", len(egfp))

    # -------------------------
    # HISTOGRAM
    # -------------------------
    plt.figure()
    plt.hist(egfp, bins=50)
    plt.title("EGFP expression (FULL DATA)")
    plt.savefig(os.path.join(outdir, f"{args.prefix}_EGFP_hist.png"), dpi=150)
    plt.close()

    # -------------------------
    # VIOLIN
    # -------------------------
    plt.figure()
    plt.violinplot(egfp, showmeans=True)
    plt.title("EGFP variability across cells")
    plt.ylabel("Expression")
    plt.savefig(os.path.join(outdir, f"{args.prefix}_EGFP_violin.png"), dpi=150)
    plt.close()

    # -------------------------
    # UMAP (if exists)
    # -------------------------
    if "X_umap" in adata.obsm:
        adata.obs["EGFP_full"] = egfp
        sc.pl.umap(adata, color="EGFP_full", show=False)
        plt.savefig(os.path.join(outdir, f"{args.prefix}_EGFP_umap.png"), dpi=150)
        plt.close()
    else:
        print("No UMAP found")

    print("DONE")


if __name__ == "__main__":
    main()

