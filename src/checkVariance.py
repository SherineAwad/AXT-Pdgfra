#!/usr/bin/env python3
import scanpy as sc
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    outdir = os.path.abspath("figures")
    os.makedirs(outdir, exist_ok=True)

    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # -------------------------
    # PCA (ONLY IF NOT EXISTS)
    # -------------------------
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=20)

    pc1 = adata.obsm["X_pca"][:, 0]
    pc2 = adata.obsm["X_pca"][:, 1]

    # -------------------------
    # VARIANCE EXPLAINED
    # -------------------------
    if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        var_ratio = adata.uns["pca"]["variance_ratio"]
    else:
        var_ratio = None

    # -------------------------
    # EGFP expression
    # -------------------------
    if "EGFP" in adata.var_names:

        egfp = adata[:, "EGFP"].X
        if hasattr(egfp, "toarray"):
            egfp = egfp.toarray().ravel()
        else:
            egfp = np.ravel(egfp)

        corr1 = np.corrcoef(pc1, egfp)[0, 1]
        corr2 = np.corrcoef(pc2, egfp)[0, 1]

        print("EGFP-PC1 correlation:", corr1)
        print("EGFP-PC2 correlation:", corr2)

        # -------------------------
        # PLOT - Boxplots
        # -------------------------
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        egfp_binary = ["EGFP+" if x > 0 else "EGFP-" for x in egfp]

        # PC1 boxplot
        pc1_neg = [pc1[i] for i in range(len(egfp_binary)) if egfp_binary[i] == "EGFP-"]
        pc1_pos = [pc1[i] for i in range(len(egfp_binary)) if egfp_binary[i] == "EGFP+"]
        bp1 = ax[0].boxplot([pc1_neg, pc1_pos], labels=["EGFP-", "EGFP+"], patch_artist=True)
        bp1['boxes'][0].set_facecolor('lightblue')
        bp1['boxes'][1].set_facecolor('lightcoral')
        ax[0].set_ylabel("PC1")
        ax[0].set_title(f"PC1 by EGFP (r={corr1:.3f})")
        ax[0].axhline(0, color='gray', linestyle='--', alpha=0.5)

        # PC2 boxplot
        pc2_neg = [pc2[i] for i in range(len(egfp_binary)) if egfp_binary[i] == "EGFP-"]
        pc2_pos = [pc2[i] for i in range(len(egfp_binary)) if egfp_binary[i] == "EGFP+"]
        bp2 = ax[1].boxplot([pc2_neg, pc2_pos], labels=["EGFP-", "EGFP+"], patch_artist=True)
        bp2['boxes'][0].set_facecolor('lightblue')
        bp2['boxes'][1].set_facecolor('lightcoral')
        ax[1].set_ylabel("PC2")
        ax[1].set_title(f"PC2 by EGFP (r={corr2:.3f})")
        ax[1].axhline(0, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()

        path = os.path.join(outdir, f"{args.prefix}_egfp_pca.png")
        plt.savefig(path, dpi=120, bbox_inches="tight")
        plt.close()

        print("Saved:", path)

    else:
        print("EGFP not found in dataset")

    # -------------------------
    # SAVE VARIANCE (IF AVAILABLE)
    # -------------------------
    if var_ratio is not None:
        df = pd.DataFrame({
            "PC": np.arange(1, len(var_ratio) + 1),
            "variance_ratio": var_ratio
        })

        csv_path = os.path.join(outdir, f"{args.prefix}_pca_variance.csv")
        df.to_csv(csv_path, index=False)

        plt.figure()
        plt.plot(df["PC"], df["variance_ratio"], marker="o")
        plt.xlabel("PC")
        plt.ylabel("Variance explained")
        plt.title("PCA variance")
        plt.savefig(os.path.join(outdir, f"{args.prefix}_variance.png"),
                    dpi=120, bbox_inches="tight")
        plt.close()

        print("Saved PCA variance outputs")

    print("DONE")

if __name__ == "__main__":
    main()
