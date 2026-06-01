import scanpy as sc
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--HVG', type=int, default=2000, 
                        help="Number of highly variable genes (0 = disable HVG filtering)")
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
    # BASIC CLEANUP
    # -------------------------
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # -------------------------
    # NORMALISATION
    # -------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()

    # -------------------------
    # HVG (OPTIONAL)
    # -------------------------
    if args.HVG > 0:
        print(f"Selecting top {args.HVG} HVGs")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.HVG,
            flavor="seurat"
        )
        adata = adata[:, adata.var.highly_variable].copy()
    else:
        print("Skipping HVG selection (using all genes)")

    # -------------------------
    # PCA
    # -------------------------
    sc.tl.pca(
        adata,
        n_comps=30,
        svd_solver="arpack"
    )

    print("PCA shape:", adata.obsm["X_pca"].shape)

    # -------------------------
    # NEIGHBORS + UMAP
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
