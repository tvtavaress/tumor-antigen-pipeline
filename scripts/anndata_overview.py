#!/usr/bin/env python3

import argparse
import html
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

MAX_EXAMPLES = 5
MAX_EXTRA_CATEGORIES = 50


# ============================================================
# HELPERS
# ============================================================

def esc(value):
    """Escape text for safe HTML rendering."""
    return html.escape(str(value))


def format_number(value):
    """Format numbers using thousands separators."""
    return f"{value:,}"


def format_bytes(size):
    """Convert bytes to a human-readable size."""
    if size is None:
        return "Unavailable"

    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(size)

    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:,.2f} {unit}"
        size /= 1024

    return "Unavailable"


def json_for_html(data):
    """Serialize data safely for an HTML script element."""
    return (
        json.dumps(data, ensure_ascii=True, default=str)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def is_missing(value):
    """Identify missing scalar metadata values."""
    try:
        result = pd.isna(value)
        return bool(result) if np.ndim(result) == 0 else False
    except (TypeError, ValueError):
        return False


def clean_value(value):
    """Convert metadata values to readable strings."""
    if is_missing(value):
        return "(missing)"

    return str(value)


def get_column_type(series):
    """Describe the pandas data type of a metadata column."""
    dtype = str(series.dtype)

    if pd.api.types.is_bool_dtype(series):
        return f"Boolean ({dtype})"

    if pd.api.types.is_integer_dtype(series):
        return f"Integer ({dtype})"

    if pd.api.types.is_float_dtype(series):
        return f"Numeric ({dtype})"

    if pd.api.types.is_datetime64_any_dtype(series):
        return f"Datetime ({dtype})"

    if isinstance(series.dtype, pd.CategoricalDtype):
        return f"Categorical ({dtype})"

    return f"Text/object ({dtype})"


def classify_column(column_name):
    """
    Infer the likely role of a metadata column from its name.
    This is a heuristic, not a guarantee of biological meaning.
    """
    name = str(column_name).lower()

    normalized = (
        name.replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )

    rules = [
        (
            "Cell type",
            [
                "cell_type",
                "celltype",
                "cell_annotation",
                "cell_ontology",
                "cell_identity",
                "cell_label",
                "cell_class",
                "cell_subtype",
                "cell_population",
                "cell_lineage",
                "annotation",
            ],
        ),
        (
            "Tumor/cancer",
            [
                "tumor_type",
                "tumour_type",
                "cancer_type",
                "cancer",
                "tumor",
                "tumour",
                "malignancy",
                "disease",
                "diagnosis",
                "histology",
                "cancer_site",
                "tumor_site",
                "tissue_type",
            ],
        ),
        (
            "Sample",
            [
                "sample",
                "specimen",
                "library",
                "sample_id",
                "sample_name",
                "sampleid",
                "barcode",
            ],
        ),
        (
            "Patient/donor",
            [
                "patient",
                "donor",
                "subject",
                "individual",
                "participant",
            ],
        ),
        (
            "Condition/treatment",
            [
                "condition",
                "treatment",
                "therapy",
                "response",
                "outcome",
                "status",
                "clinical",
                "drug",
                "exposure",
            ],
        ),
        (
            "Quality control",
            [
                "n_genes",
                "n_counts",
                "total_counts",
                "pct_counts",
                "percent",
                "mitochond",
                "doublet",
                "qc",
            ],
        ),
        (
            "Cluster",
            [
                "cluster",
                "leiden",
                "louvain",
                "resolution",
            ],
        ),
        (
            "Species",
            [
                "species",
                "organism",
                "genus",
            ],
        ),
    ]

    for category, keywords in rules:
        if any(keyword in normalized for keyword in keywords):
            return category

    return "Other metadata"


def find_candidate_columns(columns, category):
    """Return metadata columns matching an inferred category."""
    return [
        column
        for column in columns
        if classify_column(column) == category
    ]


def distribution(series):
    """Count all observed values, including missing values."""
    counts = Counter(clean_value(value) for value in series)

    return sorted(
        counts.items(),
        key=lambda item: (-item[1], item[0].lower()),
    )


def make_bar_table(counts, total, label_name="Value"):
    """Render a distribution table with proportional gray bars."""
    if not counts:
        return "<p class='muted'>No values available.</p>"

    rows = []

    for label, count in counts:
        percentage = 100 * count / total if total else 0

        rows.append(
            "<tr>"
            f"<td>{esc(label)}</td>"
            f"<td class='number'>{format_number(count)}</td>"
            f"<td class='number'>{percentage:.2f}%</td>"
            "<td class='bar-cell'>"
            "<div class='bar-track'>"
            f"<div class='bar' style='width:{min(percentage, 100):.4f}%'></div>"
            "</div>"
            "</td>"
            "</tr>"
        )

    return (
        "<div class='table-wrap'>"
        "<table class='data-table searchable'>"
        "<thead><tr>"
        f"<th>{esc(label_name)}</th>"
        "<th>Cells/records</th>"
        "<th>Percentage</th>"
        "<th>Distribution</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody>"
        "</table></div>"
    )


def make_metadata_table(frame, section_name, table_id=None):
    """Describe every column in an AnnData metadata dataframe."""
    rows = []

    for column in frame.columns:
        series = frame[column]
        missing = int(series.isna().sum())
        unique = int(series.nunique(dropna=True))

        examples = [
            clean_value(value)
            for value in series.dropna().head(MAX_EXAMPLES).tolist()
        ]

        examples_text = (
            ", ".join(examples)
            if examples
            else "No non-missing values"
        )

        rows.append(
            "<tr>"
            f"<td><code>{esc(column)}</code></td>"
            f"<td>{esc(get_column_type(series))}</td>"
            f"<td>{esc(classify_column(column))}</td>"
            f"<td class='number'>{format_number(unique)}</td>"
            f"<td class='number'>{format_number(missing)}</td>"
            f"<td class='number'>"
            f"{(100 * missing / len(frame)) if len(frame) else 0:.2f}%"
            "</td>"
            f"<td class='examples'>{esc(examples_text)}</td>"
            "</tr>"
        )

    if not rows:
        return "<p class='muted'>No columns available.</p>"

    table_id_attribute = (
        f" id='{esc(table_id)}'"
        if table_id
        else ""
    )

    return (
        f"<p class='muted'>{esc(section_name)}: "
        f"{format_number(len(frame.columns))} columns</p>"
        "<div class='table-wrap'>"
        f"<table class='data-table searchable'{table_id_attribute}>"
        "<thead><tr>"
        "<th>Column name</th>"
        "<th>Data type</th>"
        "<th>Inferred role</th>"
        "<th>Unique values</th>"
        "<th>Missing</th>"
        "<th>Missing (%)</th>"
        "<th>Example values</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody>"
        "</table></div>"
    )


# ============================================================
# HTML REPORT
# ============================================================

def build_report(adata, input_path):
    """Build a complete HTML report for an AnnData object."""

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    obs = adata.obs
    var = adata.var

    file_size = input_path.stat().st_size

    generated_at = datetime.now().astimezone().strftime(
        "%Y-%m-%d %H:%M:%S %Z"
    )

    # --------------------------------------------------------
    # Infer relevant metadata columns
    # --------------------------------------------------------

    cell_type_columns = find_candidate_columns(
        obs.columns, "Cell type"
    )

    tumor_columns = find_candidate_columns(
        obs.columns, "Tumor/cancer"
    )

    sample_columns = find_candidate_columns(
        obs.columns, "Sample"
    )

    patient_columns = find_candidate_columns(
        obs.columns, "Patient/donor"
    )

    condition_columns = find_candidate_columns(
        obs.columns, "Condition/treatment"
    )

    # --------------------------------------------------------
    # Overview metric cards
    # --------------------------------------------------------

    cards = [
        ("Cells", format_number(n_cells)),
        ("Genes", format_number(n_genes)),
        ("Cell metadata columns", format_number(len(obs.columns))),
        ("Gene metadata columns", format_number(len(var.columns))),
        ("Layers", format_number(len(adata.layers))),
        ("File size", format_bytes(file_size)),
    ]

    card_html = "".join(
        "<div class='metric-card'>"
        f"<div class='metric-label'>{esc(label)}</div>"
        f"<div class='metric-value'>{esc(value)}</div>"
        "</div>"
        for label, value in cards
    )

    # --------------------------------------------------------
    # Cell-type distributions
    # --------------------------------------------------------

    cell_type_sections = []

    for column in cell_type_columns:
        counts = distribution(obs[column])

        cell_type_sections.append(
            f"<h3>{esc(column)}</h3>"
            "<p class='muted'>"
            f"{format_number(len(counts))} distinct values, "
            "including missing values"
            "</p>"
            + make_bar_table(counts, n_cells, column)
        )

    if not cell_type_sections:
        cell_type_sections.append(
            "<p class='notice'>"
            "No cell-type column was identified automatically. "
            "Check the metadata schema below and select the correct column."
            "</p>"
        )

    # --------------------------------------------------------
    # Tumor and cancer distributions
    # --------------------------------------------------------

    tumor_sections = []

    for column in tumor_columns:
        counts = distribution(obs[column])

        tumor_sections.append(
            f"<h3>{esc(column)}</h3>"
            "<p class='muted'>"
            f"{format_number(len(counts))} distinct values, "
            "including missing values"
            "</p>"
            + make_bar_table(counts, n_cells, column)
        )

    if not tumor_sections:
        tumor_sections.append(
            "<p class='notice'>"
            "No tumor/cancer column was identified automatically. "
            "Inspect the metadata schema for the actual column name."
            "</p>"
        )

    # --------------------------------------------------------
    # Cell type x tumor cross-tabulation
    # --------------------------------------------------------

    cross_sections = []

    for cell_column in cell_type_columns:
        for tumor_column in tumor_columns:

            cross = pd.crosstab(
                obs[tumor_column].astype("string").fillna("(missing)"),
                obs[cell_column].astype("string").fillna("(missing)"),
            )

            rows = []

            for tumor, row in cross.iterrows():
                cells = "".join(
                    f"<td class='number'>{format_number(int(value))}</td>"
                    for value in row.values
                )

                rows.append(
                    f"<tr><td>{esc(tumor)}</td>{cells}</tr>"
                )

            headers = "".join(
                f"<th>{esc(cell_type)}</th>"
                for cell_type in cross.columns
            )

            cross_sections.append(
                f"<h3>{esc(tumor_column)} × {esc(cell_column)}</h3>"
                "<p class='muted'>"
                "Rows are tumor labels; columns are cell types."
                "</p>"
                "<div class='table-wrap'>"
                "<table class='data-table searchable'>"
                "<thead><tr>"
                f"<th>{esc(tumor_column)} / {esc(cell_column)}</th>"
                f"{headers}"
                "</tr></thead>"
                f"<tbody>{''.join(rows)}</tbody>"
                "</table></div>"
            )

    if not cross_sections:
        cross_sections.append(
            "<p class='muted'>"
            "A cross-tabulation requires at least one inferred "
            "cell-type column and one inferred tumor/cancer column."
            "</p>"
        )

    # --------------------------------------------------------
    # Other categorical metadata distributions
    # --------------------------------------------------------

    extra_sections = []

    excluded_columns = set(
        cell_type_columns + tumor_columns
    )

    for column in obs.columns:

        if column in excluded_columns:
            continue

        series = obs[column]
        unique_count = int(series.nunique(dropna=True))

        # Avoid displaying high-cardinality metadata such as IDs.
        if unique_count == 0 or unique_count > MAX_EXTRA_CATEGORIES:
            continue

        if not (
            pd.api.types.is_object_dtype(series)
            or pd.api.types.is_string_dtype(series)
            or isinstance(series.dtype, pd.CategoricalDtype)
            or pd.api.types.is_bool_dtype(series)
        ):
            continue

        counts = distribution(series)

        extra_sections.append(
            "<details class='distribution-detail'>"
            f"<summary>{esc(column)} "
            f"<span class='muted'>"
            f"({format_number(unique_count)} unique values)"
            "</span></summary>"
            + make_bar_table(counts, n_cells, column)
            + "</details>"
        )

    if not extra_sections:
        extra_sections.append(
            "<p class='muted'>"
            "No additional low-cardinality categorical columns found."
            "</p>"
        )

    # --------------------------------------------------------
    # AnnData structure
    # --------------------------------------------------------

    uns_items = []

    for key, value in adata.uns.items():

        try:
            description = type(value).__name__

            if isinstance(value, dict):
                description += f" ({len(value)} keys)"

            elif isinstance(value, (list, tuple, np.ndarray)):
                description += f" (length {len(value)})"

        except Exception:
            description = type(value).__name__

        uns_items.append(
            f"<tr><td><code>{esc(key)}</code></td>"
            f"<td>{esc(description)}</td></tr>"
        )

    if not uns_items:
        uns_items.append(
            "<tr><td colspan='2'>No entries in adata.uns</td></tr>"
        )

    structure = {
        "X": getattr(adata.X, "shape", None),
        "layers": {
            key: getattr(value, "shape", None)
            for key, value in adata.layers.items()
        },
        "obsm": {
            key: getattr(value, "shape", None)
            for key, value in adata.obsm.items()
        },
        "varm": {
            key: getattr(value, "shape", None)
            for key, value in adata.varm.items()
        },
        "obsp": {
            key: getattr(value, "shape", None)
            for key, value in adata.obsp.items()
        },
        "varp": {
            key: getattr(value, "shape", None)
            for key, value in adata.varp.items()
        },
    }

    structure_rows = []

    for key, value in structure.items():
        structure_rows.append(
            f"<tr><td><code>{esc(key)}</code></td>"
            f"<td><pre>{esc(json.dumps(value, indent=2, default=str))}"
            "</pre></td></tr>"
        )

    # --------------------------------------------------------
    # Contents navigation
    # --------------------------------------------------------

    contents = [
        ("overview", "01. Dataset overview"),
        ("cell-metadata", "02. Cell metadata schema — adata.obs"),
        ("cell-types", "03. Cell-type composition"),
        ("tumor-composition", "04. Tumor and cancer composition"),
        ("tumor-cell-cross", "05. Tumor × cell-type composition"),
        ("other-metadata", "06. Other categorical metadata distributions"),
    ]

    contents_html = "".join(
        "<li>"
        f"<a href='#{esc(anchor)}'>{esc(label)}</a>"
        "</li>"
        for anchor, label in contents
    )

    # --------------------------------------------------------
    # HTML document
    # --------------------------------------------------------

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">

<title>H5AD Dataset Report — {esc(input_path.name)}</title>

<style>
:root {{
    --bg: #f3f3f3;
    --panel: #ffffff;
    --panel-alt: #f8f8f8;
    --text: #252525;
    --muted: #707070;
    --border: #dedede;
    --dark: #343434;
    --mid: #858585;
    --bar-bg: #ededed;
}}

* {{
    box-sizing: border-box;
}}

html {{
    scroll-behavior: smooth;
    scroll-padding-top: 20px;
}}

body {{
    margin: 0;
    background: var(--bg);
    color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                 Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.55;
}}

header {{
    background: var(--dark);
    color: #ffffff;
    padding: 42px 6%;
}}

header h1 {{
    margin: 0 0 10px;
    font-size: clamp(25px, 4vw, 38px);
    font-weight: 650;
    letter-spacing: -0.7px;
}}

header p {{
    margin: 5px 0;
    color: #dedede;
}}

.container {{
    width: min(1500px, 94%);
    margin: 28px auto 60px;
}}

.section {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 26px;
    margin-bottom: 22px;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.025);
}}

h2 {{
    font-size: 21px;
    margin: 0 0 18px;
    padding-bottom: 12px;
    border-bottom: 1px solid var(--border);
    letter-spacing: -0.3px;
}}

h3 {{
    font-size: 16px;
    margin: 25px 0 10px;
}}

.muted {{
    color: var(--muted);
    font-size: 13px;
}}

.metrics {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(175px, 1fr));
    gap: 12px;
    margin-bottom: 22px;
}}

.metric-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 19px;
}}

.metric-label {{
    color: var(--muted);
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.8px;
}}

.metric-value {{
    font-size: 27px;
    font-weight: 650;
    margin-top: 7px;
    overflow-wrap: anywhere;
}}

.contents-list {{
    padding-left: 24px;
    margin: 0;
}}

.contents-list li {{
    padding: 5px 0;
}}

.contents-list a {{
    color: var(--text);
    text-decoration: none;
    border-bottom: 1px solid #c5c5c5;
}}

.contents-list a:hover {{
    color: #000000;
    border-bottom-color: #222222;
}}

.back-to-contents {{
    display: inline-block;
    margin-top: 18px;
    color: var(--muted);
    font-size: 13px;
    text-decoration: none;
}}

.back-to-contents:hover {{
    text-decoration: underline;
}}

.table-wrap {{
    overflow-x: auto;
    width: 100%;
    margin: 12px 0 20px;
    border: 1px solid var(--border);
    border-radius: 8px;
}}

.data-table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 13px;
    background: #ffffff;
}}

.data-table th {{
    background: #e9e9e9;
    color: #303030;
    font-weight: 650;
    text-align: left;
    white-space: nowrap;
    position: sticky;
    top: 0;
}}

.data-table th,
.data-table td {{
    padding: 11px 13px;
    border-bottom: 1px solid #e8e8e8;
    vertical-align: top;
}}

.data-table tbody tr:nth-child(even) {{
    background: #fafafa;
}}

.data-table tbody tr:hover {{
    background: #eeeeee;
}}

.number {{
    text-align: right;
    white-space: nowrap;
    font-variant-numeric: tabular-nums;
}}

.examples {{
    min-width: 200px;
    max-width: 400px;
    overflow-wrap: anywhere;
}}

code {{
    background: #eeeeee;
    padding: 2px 5px;
    border-radius: 4px;
    font-size: 12px;
    overflow-wrap: anywhere;
}}

.bar-cell {{
    min-width: 140px;
}}

.bar-track {{
    background: var(--bar-bg);
    height: 9px;
    border-radius: 20px;
    overflow: hidden;
    min-width: 90px;
    margin-top: 5px;
}}

.bar {{
    height: 100%;
    background: var(--mid);
    border-radius: 20px;
}}

.search-box {{
    width: 100%;
    max-width: 450px;
    border: 1px solid #c9c9c9;
    border-radius: 8px;
    padding: 11px 13px;
    margin: 0 0 14px;
    font-size: 14px;
    background: white;
    color: var(--text);
}}

.search-box:focus {{
    outline: 2px solid #a0a0a0;
    outline-offset: 2px;
}}

details {{
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 12px 15px;
    margin: 10px 0;
    background: var(--panel-alt);
}}

summary {{
    cursor: pointer;
    font-weight: 600;
}}

.notice {{
    background: #f1f1f1;
    border-left: 4px solid #777777;
    padding: 13px 16px;
    border-radius: 4px;
}}

pre {{
    white-space: pre-wrap;
    overflow-wrap: anywhere;
    margin: 0;
    font-size: 12px;
}}

footer {{
    text-align: center;
    color: var(--muted);
    font-size: 12px;
    padding: 15px;
}}

@media (max-width: 700px) {{
    header {{
        padding: 30px 5%;
    }}

    .section {{
        padding: 17px;
    }}

    .metric-value {{
        font-size: 23px;
    }}
}}

@media print {{
    body {{
        background: white;
    }}

    .section, .metric-card {{
        box-shadow: none;
    }}

    .search-box {{
        display: none;
    }}

    .back-to-contents {{
        display: none;
    }}
}}
</style>
</head>

<body>

<header>
    <h1>H5AD Dataset Report</h1>
    <p>{esc(input_path.name)}</p>
    <p>Generated on {esc(generated_at)}</p>
    <p>Automated structural and metadata inspection</p>
</header>

<main class="container">

    <!-- ====================================================
         CONTENTS
         Appears before Dataset overview.
         Every entry links directly to its section.
         ==================================================== -->

    <section class="section" id="contents">
        <h2>Contents</h2>
        <p>
            Select a section to navigate directly to the corresponding
            part of the report.
        </p>
        <ol class="contents-list">
            {contents_html}
        </ol>
    </section>

    <!-- ====================================================
         01. DATASET OVERVIEW
         Includes the original sections 09 and 10:
         AnnData structure and Interpretation notes.
         ==================================================== -->

    <section class="metrics">
        {card_html}
    </section>

    <section class="section" id="overview">
        <h2>01. Dataset overview</h2>

        <table class="data-table">
            <tbody>
                <tr>
                    <th>Input file</th>
                    <td>{esc(input_path.resolve())}</td>
                </tr>
                <tr>
                    <th>File size</th>
                    <td>{esc(format_bytes(file_size))}</td>
                </tr>
                <tr>
                    <th>Number of cells (observations)</th>
                    <td>{format_number(n_cells)}</td>
                </tr>
                <tr>
                    <th>Number of genes (variables)</th>
                    <td>{format_number(n_genes)}</td>
                </tr>
                <tr>
                    <th>Matrix dimensions</th>
                    <td>{esc(str(adata.shape))}</td>
                </tr>
                <tr>
                    <th>Observation index name</th>
                    <td>{esc(obs.index.name)}</td>
                </tr>
                <tr>
                    <th>Variable index name</th>
                    <td>{esc(var.index.name)}</td>
                </tr>
                <tr>
                    <th>Cell metadata columns</th>
                    <td>{format_number(len(obs.columns))}</td>
                </tr>
                <tr>
                    <th>Gene metadata columns</th>
                    <td>{format_number(len(var.columns))}</td>
                </tr>
                <tr>
                    <th>Layers</th>
                    <td>{esc(", ".join(adata.layers.keys()) or "None")}</td>
                </tr>
                <tr>
                    <th>obsm entries</th>
                    <td>{esc(", ".join(adata.obsm.keys()) or "None")}</td>
                </tr>
                <tr>
                    <th>uns entries</th>
                    <td>{esc(", ".join(adata.uns.keys()) or "None")}</td>
                </tr>
            </tbody>
        </table>

        <!-- Original section 09 moved into Dataset overview -->

        <h3 id="anndata-structure">AnnData structure</h3>

        <p>
            Matrix and embedding dimensions describe the stored structure;
            they do not guarantee that every matrix contains the same type
            of data.
        </p>

        <div class="table-wrap">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Component</th>
                        <th>Keys and dimensions</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(structure_rows)}
                </tbody>
            </table>
        </div>

        <h3>Unstructured metadata — adata.uns</h3>

        <div class="table-wrap">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Key</th>
                        <th>Stored object</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(uns_items)}
                </tbody>
            </table>
        </div>

        <!-- Original section 10 moved into Dataset overview -->

        <h3 id="interpretation-notes">Interpretation notes</h3>

        <ul>
            <li>
                Cell and tumor categories are inferred from metadata
                column names; inspect the schema to confirm the correct columns.
            </li>
            <li>
                Counts represent observations stored in this H5AD,
                not necessarily unique patients or biological replicates.
            </li>
            <li>
                Cell-type and tumor percentages use the total number
                of observations as the denominator.
            </li>
            <li>
                Missing values are counted separately and are not
                silently discarded from the main distributions.
            </li>
            <li>
                Gene counts refer to variables stored in the AnnData
                object; they do not necessarily equal the number of genes
                with nonzero expression.
            </li>
            <li>
                The report inspects metadata and object structure.
                It does not assess normalization, differential expression,
                or biological validity.
            </li>
        </ul>

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

    <!-- ====================================================
         02. CELL METADATA SCHEMA
         ==================================================== -->

    <section class="section" id="cell-metadata">
        <h2>02. Cell metadata schema — adata.obs</h2>

        <p>
            Complete inventory of observation metadata. The inferred role
            is based on column names and is intended as a guide, not a
            biological validation.
        </p>

        <input
            class="search-box"
            type="search"
            placeholder="Search metadata columns, types or examples..."
            data-table-search="obs-table"
            aria-label="Search cell metadata"
        >

        {make_metadata_table(obs, "Cell metadata", "obs-table")}

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

    <!-- ====================================================
         03. CELL-TYPE COMPOSITION
         ==================================================== -->

    <section class="section" id="cell-types">
        <h2>03. Cell-type composition</h2>

        <p>
            Counts and percentages calculated across all observations
            in the H5AD. Missing labels are displayed separately.
        </p>

        {''.join(cell_type_sections)}

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

    <!-- ====================================================
         04. TUMOR AND CANCER COMPOSITION
         ==================================================== -->

    <section class="section" id="tumor-composition">
        <h2>04. Tumor and cancer composition</h2>

        <p>
            Tumor-related metadata columns were inferred automatically
            from their names. Values are reported as stored in the dataset.
        </p>

        {''.join(tumor_sections)}

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

    <!-- ====================================================
         05. TUMOR x CELL-TYPE COMPOSITION
         ==================================================== -->

    <section class="section" id="tumor-cell-cross">
        <h2>05. Tumor × cell-type composition</h2>

        <p>
            Cross-tabulations show the number of cells for each combination
            of inferred tumor label and cell-type label.
        </p>

        {''.join(cross_sections)}

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

    <!-- ====================================================
         06. OTHER CATEGORICAL METADATA DISTRIBUTIONS
         ==================================================== -->

    <section class="section" id="other-metadata">
        <h2>06. Other categorical metadata distributions</h2>

        <p>
            Additional categorical or Boolean columns with up to
            {MAX_EXTRA_CATEGORIES} distinct non-missing values are shown
            below. All columns remain documented in the cell metadata
            schema above.
        </p>

        {''.join(extra_sections)}

        <a class="back-to-contents" href="#contents">
            ↑ Back to contents
        </a>
    </section>

</main>

<footer>
    H5AD Dataset Report · Generated with Python and AnnData · Grayscale edition
</footer>

<script>
document.querySelectorAll("[data-table-search]").forEach(function(input) {{
    input.addEventListener("input", function() {{
        const tableId = input.getAttribute("data-table-search");
        const table = document.getElementById(tableId);

        if (!table) return;

        const query = input.value.toLowerCase();

        table.querySelectorAll("tbody tr").forEach(function(row) {{
            row.style.display =
                row.textContent.toLowerCase().includes(query)
                    ? ""
                    : "none";
        }});
    }});
}});
</script>

</body>
</html>
"""


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Validate an H5AD file and generate a detailed HTML "
            "dataset and metadata report."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input H5AD file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Path for the validated output H5AD file.",
    )

    parser.add_argument(
        "--report",
        required=True,
        help="Path for the generated HTML report.",
    )

    args = parser.parse_args()

    input_path = Path(args.input).expanduser()
    output_path = Path(args.output).expanduser()
    report_path = Path(args.report).expanduser()

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input file not found: {input_path}"
        )

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Input and output paths must be different."
        )

    if report_path.resolve() in {
        input_path.resolve(),
        output_path.resolve(),
    }:
        raise ValueError(
            "The HTML report path must be different from the H5AD paths."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading H5AD: {input_path}")

    adata = ad.read_h5ad(input_path)

    print(f"Cells: {format_number(adata.n_obs)}")
    print(f"Genes: {format_number(adata.n_vars)}")
    print("Generating HTML report...")

    html_report = build_report(adata, input_path)

    report_path.write_text(
        html_report,
        encoding="utf-8",
    )

    print(f"Writing validated H5AD: {output_path}")

    adata.write_h5ad(output_path)

    print()
    print("Dataset inspection completed successfully.")
    print(f"HTML report: {report_path.resolve()}")
    print(f"Output H5AD: {output_path.resolve()}")


if __name__ == "__main__":
    main()
