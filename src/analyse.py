import scanpy as sc
import argparse
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
    adata = sc.read(args.input)

    # -------------------------
    # SAVE RAW COUNTS
    # -------------------------
    adata.layers["counts"] = adata.X.copy()

    # -------------------------
    # DOUBLETS ALREADY REMOVED (NO ACTION)
    # -------------------------
    # intentionally skipped

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # -------------------------
    # NORMALISATION
    # -------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()

    # -------------------------
    # CLEAN MATRIX
    # -------------------------
    adata.X = adata.X.copy()

    # -------------------------
    # SCALE + PCA
    # -------------------------
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(
        adata,
        n_comps=30,
        svd_solver="arpack"
    )

    print("PCA shape:", adata.obsm["X_pca"].shape)

    # -------------------------
    # UMAP
    # -------------------------
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)

    # -------------------------
    # GLOBAL UMAP
    # -------------------------
    sc.pl.umap(
        adata,
        color="sample",
        size=20,
        save=f"_{args.prefix}_umap.png"
    )

    # -------------------------
    # PER-SAMPLE UMAP
    # -------------------------
    for s in adata.obs["sample"].unique():
        sc.pl.umap(
            adata[adata.obs["sample"] == s],
            color="sample",
            size=20,
            title=f"Sample: {s}",
            save=f"_{args.prefix}_{s}.png"
        )

    # -------------------------
    # SAVE FINAL OBJECT
    # -------------------------
    adata.write(args.output)

if __name__ == "__main__":
    main()
