#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd


HPA_BASE_URL = "https://v24.proteinatlas.org/download/tsv/"

HPA_FILES = {
    "single_cell": "rna_single_cell_type.tsv.zip",
    "single_cell_tissue": "rna_single_cell_type_tissue.tsv.zip",
    "brain": "rna_single_nuclei_cluster_type.tsv.zip",
    "immune_cell": "rna_immune_cell.tsv.zip",
}

CELLTYPE_KEYWORDS = {
    "microglial cell": "central nervous system macrophage",
    "myofibroblast cell": "myofibroblast",
    "mpps": "hematopoietic",
    "cd4 t": "cd4",
    "b-cell": "mature b-cell",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Associate non-malignant candidate CTAs "
                    "with expression data from the Human Protein Atlas."
    )

    parser.add_argument(
        "--cta",
        required=True,
        help="CSV with classified candidate antigens.",
    )

    parser.add_argument(
        "--output",
        default="association.tsv",
        help="Output TSV file.",
    )

    return parser.parse_args()


def download_hpa_file(filename, directory):
    url = HPA_BASE_URL + filename
    destination = Path(directory) / filename

    print(f"Downloading: {url}")

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "tumor-antigen-pipeline/1.0"},
    )

    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with open(destination, "wb") as output:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)

    except Exception as exc:
        raise RuntimeError(
            f"Failed to download {url}: {exc}"
        ) from exc

    if not zipfile.is_zipfile(destination):
        raise RuntimeError(
            f"Downloaded file is not a valid ZIP archive: {url}"
        )

    print(f"Downloaded: {filename}")

    return destination


def read_hpa_zip(path):
    with zipfile.ZipFile(path) as archive:
        members = [
            name for name in archive.namelist()
            if not name.endswith("/")
            and name.lower().endswith((".tsv", ".txt"))
        ]

        if len(members) != 1:
            raise ValueError(
                f"Expected one TSV/TXT file in {path}; "
                f"found: {members}"
            )

        with archive.open(members[0]) as handle:
            return pd.read_csv(handle, sep="\t", low_memory=False)


def normalize_gene(series):
    return series.astype("string").str.strip().str.upper()


def load_cta(path):
    cta = pd.read_csv(path)

    required = {
        "type",
        "specificity_profile",
        "celltypist_labels",
        "gene_name",
        "pct_expressing",
    }

    missing = required - set(cta.columns)

    if missing:
        raise ValueError(
            f"Missing columns in CTA CSV: {sorted(missing)}"
        )

    cta = cta[
        cta["type"].astype("string").str.strip().str.lower()
        == "cancer-testis"
    ].copy()

    cta = cta[
        ~cta["specificity_profile"]
        .astype("string")
        .str.strip()
        .str.lower()
        .isin(["low specificity", "not detected"])
    ].copy()

    cta = cta[
        cta["celltypist_labels"]
        .astype("string")
        .str.strip()
        .str.lower()
        != "malignant cell"
    ].copy()

    cta["gene_name"] = normalize_gene(cta["gene_name"])

    cta = cta[
        cta["gene_name"].notna()
        & cta["gene_name"].ne("")
        & cta["gene_name"].ne("NAN")
    ].copy()

    return cta


def expression_column(df, candidates, dataset_name):
    for column in candidates:
        if column in df.columns:
            return column

    raise ValueError(
        f"No supported expression column found in {dataset_name}. "
        f"Expected one of {candidates}; "
        f"available columns: {list(df.columns)}"
    )


def require_columns(df, columns, dataset_name):
    missing = set(columns) - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns in {dataset_name}: {sorted(missing)}"
        )


def make_distribution(datasets, genes):
    distribution = []

    # 1. Single-cell type
    df = datasets["single_cell"].copy()
    require_columns(df, ["Gene name", "Cell type"], "HPA single cell")

    metric_col = expression_column(
        df, ["NX", "nTPM"], "HPA single cell"
    )

    df["Gene name"] = normalize_gene(df["Gene name"])
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df[df["Gene name"].isin(genes)]

    distribution.append(pd.DataFrame({
        "Gene name": df["Gene name"],
        "HPA dataset": "Single Cell",
        "context": df["Cell type"].astype("string"),
        "metric": metric_col,
        "expression": df[metric_col],
    }))

    # 2. Single-cell tissue
    df = datasets["single_cell_tissue"].copy()
    require_columns(
        df, ["Gene name", "Tissue", "Cell type"],
        "HPA single cell tissue"
    )

    metric_col = expression_column(
        df, ["pTPM", "nTPM"], "HPA single cell tissue"
    )

    df["Gene name"] = normalize_gene(df["Gene name"])
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df[df["Gene name"].isin(genes)]

    distribution.append(pd.DataFrame({
        "Gene name": df["Gene name"],
        "HPA dataset": "Single Cell Tissue",
        "context": (
            df["Tissue"].astype("string")
            + " | "
            + df["Cell type"].astype("string")
        ),
        "metric": metric_col,
        "expression": df[metric_col],
    }))

    # 3. Single-nuclei brain
    df = datasets["brain"].copy()
    require_columns(
        df, ["Gene name", "Cluster type"], "HPA brain"
    )

    metric_col = expression_column(
        df, ["nTPM", "NX"], "HPA brain"
    )

    df["Gene name"] = normalize_gene(df["Gene name"])
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df[df["Gene name"].isin(genes)]

    distribution.append(pd.DataFrame({
        "Gene name": df["Gene name"],
        "HPA dataset": "Single Nuclei Brain",
        "context": df["Cluster type"].astype("string"),
        "metric": metric_col,
        "expression": df[metric_col],
    }))

    # 4. Immune cells
    df = datasets["immune_cell"].copy()
    require_columns(
        df, ["Gene name", "Immune cell"], "HPA immune cell"
    )

    metric_col = expression_column(
        df, ["pTPM", "nTPM"], "HPA immune cell"
    )

    df["Gene name"] = normalize_gene(df["Gene name"])
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df[df["Gene name"].isin(genes)]

    distribution.append(pd.DataFrame({
        "Gene name": df["Gene name"],
        "HPA dataset": "Immune Cell",
        "context": df["Immune cell"].astype("string").str.strip(),
        "metric": metric_col,
        "expression": df[metric_col],
    }))

    return pd.concat(distribution, ignore_index=True)


def associate_ctas(cta, distribution):
    associations = []

    distribution["context"] = (
        distribution["context"].fillna("").astype(str)
    )

    for _, row in cta.iterrows():
        gene = row["gene_name"]

        celltype = str(row["celltypist_labels"]).strip().lower()

        if not celltype or celltype == "nan":
            continue

        subset = distribution[
            distribution["Gene name"] == gene
        ].copy()

        if subset.empty:
            continue

        pattern = rf"\b{re.escape(celltype)}\b"

        matched = subset[
            subset["context"].str.contains(
                pattern, case=False, na=False, regex=True
            )
        ].copy()

        if matched.empty and celltype in CELLTYPE_KEYWORDS:
            keyword = CELLTYPE_KEYWORDS[celltype]

            matched = subset[
                subset["context"].str.contains(
                    re.escape(keyword), case=False, na=False
                )
            ].copy()

        if matched.empty:
            continue

        matched["celltypist_labels"] = row["celltypist_labels"]
        matched["pct_expressing"] = row["pct_expressing"]
        matched["specificity_profile"] = row["specificity_profile"]

        associations.append(matched)

    columns = [
        "Gene name",
        "HPA dataset",
        "context",
        "metric",
        "expression",
        "celltypist_labels",
        "pct_expressing",
        "specificity_profile",
    ]

    if not associations:
        return pd.DataFrame(columns=columns)

    return pd.concat(
        associations, ignore_index=True
    )[columns]


def main():
    args = parse_args()

    print("Loading candidate CTAs...")
    cta = load_cta(args.cta)
    genes = set(cta["gene_name"].dropna().unique())

    print(f"Unique candidate genes: {len(genes)}")

    # Temporary directory is automatically deleted after processing.
    with tempfile.TemporaryDirectory(prefix="hpa_download_") as tmpdir:
        downloaded = {}

        for key, filename in HPA_FILES.items():
            zip_path = download_hpa_file(filename, tmpdir)
            downloaded[key] = read_hpa_zip(zip_path)

        print("Building expression associations...")
        distribution = make_distribution(downloaded, genes)
        association = associate_ctas(cta, distribution)

        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)

        association.to_csv(
            output,
            sep="\t",
            index=False,
        )

    print(f"Associations generated: {len(association)}")
    print(f"Output: {args.output}")
    print("Temporary HPA downloads have been removed.")


if __name__ == "__main__":
    main()
