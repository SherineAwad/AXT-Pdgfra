import scanpy as sc
import argparse
import os
import harmonypy as hm
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    # -------------------------
    # LOAD DATA
    # -------------------------
    adata = sc.read_h5ad(args.input)

    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # -------------------------
    # BATCH KEY
    # -------------------------
    if "batch" in adata.obs.columns:
        batch_key = "batch"
    elif "sample" in adata.obs.columns:
        batch_key = "sample"
    else:
        raise ValueError("No batch column found")

    print(f"Using batch key: {batch_key}")

    # -------------------------
    # HARMONY FIX
    # -------------------------
    X_pca = adata.obsm["X_pca"]

    meta_data = adata.obs[batch_key].values

    # REQUIRED FIX: vars_use (dummy is fine)
    vars_use = np.arange(X_pca.shape[1])

    ho = hm.run_harmony(X_pca, meta_data, vars_use)

    # FIX: Correct the transpose - ho.Z_corr already has correct orientation (cells × PCs)
    adata.obsm["X_pca_harmony"] = ho.Z_corr

    # -------------------------
    # UMAP
    # -------------------------
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)

    # -------------------------
    # PLOTS
    # -------------------------
    os.makedirs("figures", exist_ok=True)

    sc.pl.umap(
        adata,
        color=batch_key,
        save=f"_{args.prefix}_harmony_batch.png"
    )

    sc.pl.umap(
        adata,
        color="sample" if "sample" in adata.obs.columns else batch_key,
        save=f"_{args.prefix}_harmony_sample.png"
    )

    # -------------------------
    # SAVE
    # -------------------------
    adata.write(args.output, compression="gzip")

    print(f"SAVED → {args.output}")

if __name__ == "__main__":
    main()
