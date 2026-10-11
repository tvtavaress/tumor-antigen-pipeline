
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def calculate_expression_metrics(
    input_h5ad,
    celltype_col="cell_type_grouped",
    gene_col="feature_name"
):
    """Calculate expression metrics without saving an intermediate CSV."""

    adata = sc.read_h5ad(input_h5ad)

    if celltype_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{celltype_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    if gene_col in adata.var.columns:
        gene_names = adata.var[gene_col].astype(str).to_numpy()
    else:
        gene_names = adata.var_names.astype(str).to_numpy()

    cell_types = adata.obs[celltype_col]
    results = []

    for cell_type in cell_types.dropna().unique():
        mask = (cell_types == cell_type).to_numpy()
        matrix = adata.X[mask, :]
        n_cells = matrix.shape[0]

        if n_cells == 0:
            continue

        if sp.issparse(matrix):
            n_expressing = np.asarray(
                (matrix > 0).sum(axis=0)
            ).ravel()
        else:
            n_expressing = np.count_nonzero(
                np.asarray(matrix) > 0,
                axis=0
            )

        pct_expressing = (n_expressing / n_cells) * 100

        results.append(
            pd.DataFrame({
                "gene_name": gene_names,
                "celltypist_labels": str(cell_type),
                "pct_expressing": pct_expressing
            })
        )

    if not results:
        raise ValueError(
            f"No cell types found in '{celltype_col}'."
        )

    return pd.concat(results, ignore_index=True)


def classify_gene_fraction(
    subdf,
    min_detect,
    group_frac,
    fold
):
    eps = 1e-12

    subdf = subdf.sort_values(
        "pct_expressing",
        ascending=False
    ).copy()

    top = subdf.iloc[0]["pct_expressing"]
    top_celltype = subdf.iloc[0]["celltypist_labels"]

    second = (
        subdf.iloc[1]["pct_expressing"]
        if len(subdf) > 1
        else 0.0
    )

    others = subdf.iloc[1:]["pct_expressing"]
    others_mean = others.mean() if len(others) > 0 else 0.0

    if top < min_detect:
        subdf["specificity_profile"] = "Not detected"
        return subdf

    subdf["specificity_profile"] = "Low specificity"

    subdf.loc[
        subdf["pct_expressing"] < min_detect,
        "specificity_profile"
    ] = "Not detected"

    if top >= fold * max(second, eps):
        subdf.loc[
            subdf["celltypist_labels"] == top_celltype,
            "specificity_profile"
        ] = "Enriched"
        return subdf

    high_group = subdf["pct_expressing"] >= group_frac * top

    if 2 <= high_group.sum() <= 5:
        mean_group = subdf.loc[
            high_group, "pct_expressing"
        ].mean()

        mean_rest = subdf.loc[
            ~high_group, "pct_expressing"
        ].mean()

        if mean_group >= fold * mean_rest:
            subdf.loc[
                high_group,
                "specificity_profile"
            ] = "Group enriched"
            return subdf

    if top >= fold * others_mean:
        subdf.loc[
            subdf["celltypist_labels"] == top_celltype,
            "specificity_profile"
        ] = "Enhanced"
        return subdf

    return subdf


def classify_specificity(
    df,
    min_detect=1.0,
    group_frac=0.25,
    fold=4
):
    required_cols = {
        "gene_name",
        "celltypist_labels",
        "pct_expressing"
    }

    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    classified_groups = []

    for _, subdf in df.groupby("gene_name", sort=False):
        classified_groups.append(
            classify_gene_fraction(
                subdf,
                min_detect=min_detect,
                group_frac=group_frac,
                fold=fold
            )
        )

    return pd.concat(classified_groups, ignore_index=True)


def save_outputs(
    classified,
    output_dir,
    target_celltype="malignant cell"
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d-%H%M")

    classified_path = (
        output_dir / f"classified_specificity_{timestamp}.csv"
    )

    summary_path = (
        output_dir / f"specificity_summary_{timestamp}.csv"
    )

    target_path = (
        output_dir
        / f"enriched_{target_celltype.replace(' ', '_')}_{timestamp}.csv"
    )

    classified.to_csv(classified_path, index=False)

    summary = (
        classified.groupby(
            ["celltypist_labels", "specificity_profile"]
        )["gene_name"]
        .nunique()
        .reset_index()
        .rename(columns={"gene_name": "n_genes"})
    )

    summary.to_csv(summary_path, index=False)

    enriched_target = classified[
        (classified["celltypist_labels"] == target_celltype)
        & (classified["specificity_profile"] == "Enriched")
    ]

    enriched_target.to_csv(target_path, index=False)

    return {
        "classified_specificity": classified_path,
        "specificity_summary": summary_path,
        "enriched_target": target_path
    }


def run_specificity(
    input_h5ad,
    output_dir,
    min_detect=1.0,
    group_frac=0.25,
    fold=4,
    target_celltype="malignant cell"
):
    metrics = calculate_expression_metrics(input_h5ad)

    classified = classify_specificity(
        metrics,
        min_detect=min_detect,
        group_frac=group_frac,
        fold=fold
    )

    return save_outputs(
        classified,
        output_dir,
        target_celltype=target_celltype
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate expression metrics in memory and "
            "classify cell-type expression specificity."
        )
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)

    args = parser.parse_args()

    outputs = run_specificity(
        input_h5ad=args.input,
        output_dir=args.output_dir
    )

    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
