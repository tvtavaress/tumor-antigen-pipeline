#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate candidate-gene expression percentages across "
            "malignant cells from multiple tumor types and generate "
            "a binary tumor-expression matrix."
        )
    )

    parser.add_argument("--input", required=True, help="Input H5AD file.")
    parser.add_argument(
        "--gene-name-mapping",
        default=None,
        help="Optional JSON mapping Ensembl IDs to gene symbols.",
    )
    parser.add_argument(
        "--candidate-genes",
        required=True,
        help="CSV containing candidate genes.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/candidates_per_tumors",
    )
    parser.add_argument("--celltype-column", default="cell_type")
    parser.add_argument("--malignant-label", default="malignant cell")
    parser.add_argument("--tumor-column", default="disease")
    parser.add_argument("--candidate-gene-column", default="gene_name")
    parser.add_argument("--expression-threshold", type=float, default=0.1)
    parser.add_argument("--binary-threshold", type=float, default=1.0)

    return parser.parse_args()


def save_histogram(values, output_path, xlabel, title, bins):
    """Save a histogram, including when no values are available."""
    values = np.asarray(values, dtype=float).ravel()
    values = values[np.isfinite(values)]

    plt.figure(figsize=(7, 5))

    if values.size > 0:
        plt.hist(values, bins=bins, edgecolor="black")
    else:
        plt.text(
            0.5,
            0.5,
            "No candidate genes found",
            ha="center",
            va="center",
            transform=plt.gca().transAxes,
        )
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        print(f"WARNING: No values available for histogram: {title}")

    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Saved: {output_path}")


def load_gene_mapping(mapping_path):
    """
    Load an Ensembl ID -> gene symbol mapping and create its inverse.

    Example:
        ENSG00000106541 -> AGR2
        AGR2 -> ENSG00000106541
    """
    if not mapping_path:
        print(
            "No gene-name mapping provided. "
            "Candidate matching will use the values as supplied."
        )
        return None, None

    with open(mapping_path, "r", encoding="utf-8") as handle:
        ensembl_to_symbol = json.load(handle)

    if not isinstance(ensembl_to_symbol, dict):
        raise ValueError(
            "The gene-name mapping JSON must contain a dictionary "
            "mapping Ensembl IDs to gene symbols."
        )

    # Normalize mapping values and keys to strings.
    ensembl_to_symbol = {
        str(ensembl_id).strip(): str(symbol).strip()
        for ensembl_id, symbol in ensembl_to_symbol.items()
        if ensembl_id is not None and symbol is not None
    }

    # Invert the mapping: gene symbol -> Ensembl ID.
    # If a symbol occurs more than once, keep the first mapping.
    symbol_to_ensembl = {}

    for ensembl_id, symbol in ensembl_to_symbol.items():
        if symbol and symbol not in symbol_to_ensembl:
            symbol_to_ensembl[symbol] = ensembl_id

    print(f"Loaded gene-name mapping: {mapping_path}")
    print(f"Ensembl IDs in mapping: {len(ensembl_to_symbol):,}")
    print(f"Gene symbols in inverse mapping: {len(symbol_to_ensembl):,}")

    return ensembl_to_symbol, symbol_to_ensembl


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load Ensembl <-> gene-symbol mappings.
    ensembl_to_symbol, symbol_to_ensembl = load_gene_mapping(
        args.gene_name_mapping
    )

    # Load candidate genes from the CSV.
    candidates_df = pd.read_csv(args.candidate_genes)

    if args.candidate_gene_column not in candidates_df.columns:
        raise ValueError(
            f"Column '{args.candidate_gene_column}' not found in "
            f"{args.candidate_genes}. Available columns: "
            f"{list(candidates_df.columns)}"
        )

    candidate_genes = (
        candidates_df[args.candidate_gene_column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    candidate_genes = [
        gene for gene in candidate_genes.drop_duplicates().tolist()
        if gene
    ]

    print(f"Candidate genes in input file: {len(candidate_genes):,}")

    # Load H5AD.
    print("Loading H5AD...")
    adata = sc.read_h5ad(args.input)

    for column in (args.celltype_column, args.tumor_column):
        if column not in adata.obs.columns:
            raise ValueError(
                f"Column '{column}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    n_total = adata.n_obs

    # Filter malignant cells.
    adata = adata[
        adata.obs[args.celltype_column].astype(str)
        == args.malignant_label
    ].copy()

    n_malignant = adata.n_obs

    if n_malignant == 0:
        raise ValueError(
            f"No malignant cells found using "
            f"{args.celltype_column} == '{args.malignant_label}'."
        )

    print(f"Total cells: {n_total:,}")
    print(f"Malignant cells: {n_malignant:,}")
    print(f"Percentage: {100 * n_malignant / n_total:.2f}%")

    # Calculate expression percentages for every gene in each tumor.
    tumors = sorted(
        adata.obs[args.tumor_column].dropna().unique().tolist(),
        key=str,
    )

    genes = adata.var_names.astype(str)
    result = pd.DataFrame(index=genes)
    result.index.name = "Gene_ID"

    for tumor in tumors:
        print(f"Processing tumor: {tumor}")

        mask = (
            adata.obs[args.tumor_column].to_numpy() == tumor
        )
        sub = adata[mask]

        if sub.n_obs == 0:
            continue

        X = sub.X

        if sparse.issparse(X):
            pct = (
                np.asarray(
                    (X > args.expression_threshold).mean(axis=0)
                ).ravel()
                * 100
            )
        else:
            X = np.asarray(X)
            pct = (
                (X > args.expression_threshold).mean(axis=0) * 100
            )

        result[str(tumor)] = pct

    # Histogram for all genes.
    save_histogram(
        result.to_numpy().ravel(),
        output_dir / "histogram_all_genes.png",
        "% malignant cells expressing gene",
        "Distribution of gene detection across tumor types",
        bins=40,
    )

    # Match candidate symbols to Ensembl IDs in the H5AD.
    available_genes = set(result.index)

    candidate_gene_to_ensembl = {}

    for gene in candidate_genes:
        if gene in available_genes:
            # Also support a candidate CSV that already contains Ensembl IDs.
            candidate_gene_to_ensembl[gene] = gene
        elif symbol_to_ensembl is not None:
            ensembl_id = symbol_to_ensembl.get(gene)

            if ensembl_id is not None:
                candidate_gene_to_ensembl[gene] = ensembl_id
            else:
                candidate_gene_to_ensembl[gene] = None
        else:
            candidate_gene_to_ensembl[gene] = None

    candidate_genes_found = []
    candidate_genes_missing = []

    # Preserve the CSV order and avoid duplicate Ensembl IDs.
    seen_ensembl_ids = set()

    for gene in candidate_genes:
        ensembl_id = candidate_gene_to_ensembl[gene]

        if ensembl_id in available_genes:
            if ensembl_id not in seen_ensembl_ids:
                candidate_genes_found.append((gene, ensembl_id))
                seen_ensembl_ids.add(ensembl_id)
        else:
            candidate_genes_missing.append(gene)

    print(
        f"Candidate genes found in H5AD: "
        f"{len(candidate_genes_found):,}"
    )
    print(
        f"Candidate genes absent from H5AD: "
        f"{len(candidate_genes_missing):,}"
    )

    if candidate_genes_missing:
        print(
            "WARNING: Some candidate genes could not be matched to "
            "Ensembl IDs present in adata.var_names."
        )
        print(
            "First missing candidates: "
            + ", ".join(candidate_genes_missing[:20])
        )

    # Select the matched Ensembl IDs and retain their gene symbols.
    matched_ensembl_ids = [
        ensembl_id for _, ensembl_id in candidate_genes_found
    ]

    result_candidates = result.loc[matched_ensembl_ids].copy()
    result_candidates = result_candidates.reset_index()

    # Add gene symbols using the original Ensembl -> symbol mapping.
    if ensembl_to_symbol is not None:
        result_candidates.insert(
            1,
            "Gene_Name",
            result_candidates["Gene_ID"].map(ensembl_to_symbol),
        )
    else:
        # Without a mapping file, retain the candidate's original label
        # whenever it can be matched directly.
        ensembl_to_candidate = {
            ensembl_id: symbol
            for symbol, ensembl_id in candidate_genes_found
        }
        result_candidates.insert(
            1,
            "Gene_Name",
            result_candidates["Gene_ID"].map(ensembl_to_candidate),
        )

    # Save expression percentages.
    percentage_path = (
        output_dir / "candidate_gene_pct_expression_by_tumor.csv"
    )

    result_candidates.to_csv(percentage_path, index=False)
    print(f"Saved: {percentage_path}")

    # Identify tumor columns regardless of whether Gene_Name exists.
    tumor_columns = [
        column
        for column in result_candidates.columns
        if column not in ("Gene_ID", "Gene_Name")
    ]

    # Histogram for candidate genes.
    save_histogram(
        result_candidates[tumor_columns].to_numpy().ravel(),
        output_dir / "histogram_candidate_genes.png",
        "% malignant cells expressing candidate gene",
        "Distribution of candidate gene detection",
        bins=30,
    )

    # Generate binary expression matrix.
    binary = result_candidates.copy()

    if tumor_columns:
        binary[tumor_columns] = (
            binary[tumor_columns] >= args.binary_threshold
        ).astype(int)

    binary_path = output_dir / "candidate_gene_binary_matrix.csv"
    binary.to_csv(binary_path, index=False)

    print(f"Saved: {binary_path}")

    # Save missing candidates to facilitate debugging.
    if candidate_genes_missing:
        missing_path = output_dir / "candidate_genes_not_found.csv"

        pd.DataFrame(
            {"Gene_Name": candidate_genes_missing}
        ).to_csv(missing_path, index=False)

        print(f"Saved missing candidates list: {missing_path}")

    print("\nAnalysis completed successfully.")


if __name__ == "__main__":
    main()
