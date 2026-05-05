import scanpy as sc
import argparse
import os
import harmonypy as hm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Run Harmony as in tutorial
    pcs = adata.obsm['X_pca']
    harmony_out = hm.run_harmony(pcs, adata.obs, "batch")
    adata.obsm['X_pca_harmony'] = harmony_out.Z_corr

    # UMAP
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)

    # Plot
    os.makedirs('figures', exist_ok=True)
    sc.pl.umap(adata, color='batch', save=f'_{args.prefix}_harmony_batch.png')

    # Save
    adata.write(args.output, compression='gzip')
    print(f'SAVED → {args.output}')

if __name__ == '__main__':
    main()
