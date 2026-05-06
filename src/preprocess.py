import scanpy as sc
import sys
import importlib_metadata
import argparse
import numpy as np
import scipy.sparse as sp

sys.modules['importlib.metadata'] = importlib_metadata


def main():
    parser = argparse.ArgumentParser()

    # Required args
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)

    # QC parameters
    parser.add_argument('--min_genes', type=int, default=800)
    parser.add_argument('--max_genes', type=int, default=8000)
    parser.add_argument('--min_counts', type=int, default=1200)
    parser.add_argument('--max_counts', type=int, default=50000)
    parser.add_argument('--max_mt', type=float, default=10)

    parser.add_argument('--min_cells_gene', type=int, default=3)
    parser.add_argument('--min_genes_cell', type=int, default=100)

    args = parser.parse_args()

    # ----------------------------
    # LOAD H5AD
    # ----------------------------
    print(f"Reading {args.input}...")
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    # ----------------------------
    # QC METRICS (MT + RIBO)
    # ----------------------------

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")

    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.lower().str.startswith(("rpl", "rps", "mrpl", "mrps"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        inplace=True,
        log1p=True
    )

    # ----------------------------
    # BEFORE QC PLOTS
    # ----------------------------
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{args.prefix}_preQC.png"
    )

    # -----------------------------
    # FILTER CELLS
    # -----------------------------
    adata = adata[
        (adata.obs['n_genes_by_counts'] > args.min_genes) &
        (adata.obs['n_genes_by_counts'] < args.max_genes) &
        (adata.obs['total_counts'] > args.min_counts) &
        (adata.obs['total_counts'] < args.max_counts) &
        (adata.obs['pct_counts_mt'] < args.max_mt), :
    ]

    sc.pp.filter_cells(adata, min_genes=args.min_genes_cell)
    sc.pp.filter_genes(adata, min_cells=args.min_cells_gene)

    # -----------------------------
    # REMOVE ZERO EXPRESSION GENES
    # -----------------------------
    if sp.issparse(adata.X):
        gene_sums = np.array(adata.X.sum(axis=0)).flatten()
    else:
        gene_sums = adata.X.sum(axis=0)

    adata = adata[:, gene_sums > 0].copy()

    # ----------------------------
    # AFTER QC PLOTS
    # ----------------------------
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{args.prefix}_AfterQC.png"
    )

    # ----------------------------
    # SAVE
    # ----------------------------
    adata.write(args.output)


if __name__ == "__main__":
    main()
