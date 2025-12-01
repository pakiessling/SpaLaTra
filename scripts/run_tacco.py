import argparse

import numpy as np
import scanpy as sc
import tacco as tc


def label_transfer(adata, ref, ct_column="cell_type", res_column="tacco", **kwargs):
    assert adata.X.max().is_integer() and ref.X.max().is_integer(), "Data is transformed!"
    result_df = tc.tl.annotate(
        adata,
        ref,
        annotation_key=ct_column,
        result_key=None,
        # skip_checks=True,
        **kwargs,
    )

    return result_df


def main():
    parse = argparse.ArgumentParser("Label transfer with Tacco")
    parse.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the folder containing input h5ad files",
    )
    parse.add_argument("--ref", type=str, required=True, help="Path to the reference h5ad file")
    parse.add_argument(
        "--mc",
        type=int,
        required=False,
        help="How many multi centers to use, equal to batches in ref",
    )
    parse.add_argument(
        "--ct_column", type=str, default="cell_subtype", help="Column in ref where annotation is"
    )
    parse.add_argument(
        "--res_column", type=str, help="Column where result should be saved", default="tacco"
    )
    parse.add_argument(
        "--layer",
        type=str,
        help="What layer in the reference to use",
    )
    parse.add_argument(
        "--output",
        type=str,
        help="Where to save result csv",
    )
    args = parse.parse_args()

    # Load reference data
    print(f"Loading reference data from {args.ref}")
    ref = sc.read_h5ad(args.ref)

    if args.layer:
        ref.X = ref.layers[args.layer]

    # Load the full dataset
    print(f"Loading input data from {args.input}")
    full_adata = sc.read_h5ad(args.input)
    sc.pp.filter_genes(full_adata, min_counts=1)
    from scipy import sparse

    if sparse.issparse(full_adata.X):
        full_adata.X.data = np.round(full_adata.X.data).astype(np.float64)
    else:
        full_adata.X = np.round(full_adata.X).astype(np.float64)

    # Set multi_center parameter
    args.mc = ref.obs.Donor.nunique() if args.mc is None else args.mc

    final_results = label_transfer(
        full_adata, ref, ct_column=args.ct_column, res_column=args.res_column, multi_center=args.mc
    )
    # Save the results
    print(f"Saving results to {args.output}")
    final_results.to_csv(args.output)
    print("Processing complete")


if __name__ == "__main__":
    main()
