import argparse
import glob
import os

import anndata
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

parser = argparse.ArgumentParser(description="Generate SPaLaTra QC report")
parser.add_argument("--consensus", required=True, help="Path to consensus.csv")
parser.add_argument("--input", required=True, help="Directory of query .h5ad files, or path to a single .h5ad file")
parser.add_argument("--output", required=True, help="Output HTML path")
parser.add_argument("--embedding", default="spatial", help="Unused; kept for CLI compatibility")
args = parser.parse_args()

# ── Map cells to samples ───────────────────────────────────────────────────────
if os.path.isfile(args.input) and args.input.endswith(".h5ad"):
    h5ad_files = [args.input]
else:
    h5ad_files = glob.glob(os.path.join(args.input, "*.h5ad"))
if not h5ad_files:
    raise FileNotFoundError(f"No .h5ad files found in {args.input}")

cell_sample = {}
for path in h5ad_files:
    sample_name = os.path.splitext(os.path.basename(path))[0]
    adata = anndata.read_h5ad(path, backed="r")
    for cell in adata.obs_names:
        cell_sample[cell] = sample_name
    adata.file.close()

# ── Load consensus ─────────────────────────────────────────────────────────────
consensus = pd.read_csv(args.consensus, index_col=0)
consensus["sample"] = consensus.index.map(cell_sample)

METHOD_COLS = ["tacco", "singler", "rctd", "phispace", "insitutype", "consensus"]
PRESENT_METHODS = [c for c in METHOD_COLS if c in consensus.columns]
PRIMARY_METHODS = [c for c in ["tacco", "singler", "rctd", "phispace", "insitutype"] if c in consensus.columns]

figures = []

# ── Section 1: Pairwise method agreement heatmap ──────────────────────────────
n = len(PRIMARY_METHODS)
agreement_matrix = np.zeros((n, n))
for i, m1 in enumerate(PRIMARY_METHODS):
    for j, m2 in enumerate(PRIMARY_METHODS):
        if m1 == m2:
            agreement_matrix[i, j] = 100.0
            continue
        valid = consensus[[m1, m2]].dropna()
        valid = valid[(valid[m1] != "?") & (valid[m2] != "?")]
        if len(valid) == 0:
            agreement_matrix[i, j] = float("nan")
        else:
            agreement_matrix[i, j] = (valid[m1] == valid[m2]).mean() * 100

agreement_df = pd.DataFrame(agreement_matrix, index=PRIMARY_METHODS, columns=PRIMARY_METHODS)
fig1 = px.imshow(
    agreement_df,
    text_auto=".1f",
    color_continuous_scale="Blues",
    zmin=0,
    zmax=100,
    title="Pairwise Method Agreement (% cells with same label)",
    labels=dict(color="Agreement (%)"),
)
fig1.update_layout(height=500)
figures.append(("Pairwise Method Agreement", fig1))

# ── Section 2: Per-sample agreement ───────────────────────────────────────────
if "agreement_score" in consensus.columns and "sample" in consensus.columns:
    sample_agg = (
        consensus.groupby("sample")["agreement_score"]
        .agg(mean_agreement="mean", n_cells="count")
        .reset_index()
        .sort_values("mean_agreement")
    )
    fig2 = px.bar(
        sample_agg,
        x="sample",
        y="mean_agreement",
        hover_data=["n_cells"],
        title="Mean Agreement Score per Sample (sorted worst → best)",
        labels={"mean_agreement": "Mean Agreement Score", "sample": "Sample", "n_cells": "Cells"},
        color="mean_agreement",
        color_continuous_scale="RdYlGn",
        range_color=[0, 1],
    )
    fig2.update_layout(height=500, xaxis_tickangle=-45, coloraxis_showscale=False)
    figures.append(("Per-Sample Agreement", fig2))

# ── Section 3: Per-cell-type agreement ────────────────────────────────────────
if "agreement_score" in consensus.columns and "consensus" in consensus.columns:
    ct_agg = (
        consensus.groupby("consensus")["agreement_score"]
        .agg(mean_agreement="mean", n_cells="count")
        .reset_index()
        .sort_values("mean_agreement")
    )
    fig3 = px.bar(
        ct_agg,
        x="consensus",
        y="mean_agreement",
        hover_data=["n_cells"],
        title="Mean Agreement Score per Cell Type (sorted worst → best)",
        labels={"mean_agreement": "Mean Agreement Score", "consensus": "Cell Type", "n_cells": "Cells"},
        color="mean_agreement",
        color_continuous_scale="RdYlGn",
        range_color=[0, 1],
    )
    fig3.update_layout(height=500, xaxis_tickangle=-45, coloraxis_showscale=False)
    figures.append(("Per-Cell-Type Agreement", fig3))

# ── Section 4: Label frequency bar chart ──────────────────────────────────────
freq_rows = []
for method in PRESENT_METHODS:
    counts = consensus[method].value_counts()
    for label, count in counts.items():
        freq_rows.append({"Method": method, "Label": label, "Count": int(count)})
freq_df = pd.DataFrame(freq_rows)

fig4 = px.bar(
    freq_df,
    x="Label",
    y="Count",
    color="Method",
    barmode="group",
    title="Cell Type Label Frequency per Method",
)
fig4.update_layout(height=500, xaxis_tickangle=-45)
figures.append(("Label Frequency", fig4))

# ── Section 5: Quality metrics ────────────────────────────────────────────────
fig5 = make_subplots(rows=1, cols=2, subplot_titles=["Agreement Score Distribution", "Ambiguous Cell Count"])

if "agreement_score" in consensus.columns:
    scores = consensus["agreement_score"].dropna()
    fig5.add_trace(
        go.Histogram(x=scores, nbinsx=20, name="Agreement Score", marker_color="#636EFA"),
        row=1, col=1,
    )
    fig5.update_xaxes(title_text="Agreement Score (0–1)", row=1, col=1)
    fig5.update_yaxes(title_text="Cell Count", row=1, col=1)

if "is_ambiguous" in consensus.columns:
    ambig_counts = consensus["is_ambiguous"].value_counts().reset_index()
    ambig_counts.columns = ["Ambiguous", "Count"]
    ambig_counts["Ambiguous"] = ambig_counts["Ambiguous"].map({True: "Ambiguous", False: "Resolved"})
    fig5.add_trace(
        go.Bar(
            x=ambig_counts["Ambiguous"],
            y=ambig_counts["Count"],
            name="Ambiguity",
            marker_color=["#EF553B", "#00CC96"],
        ),
        row=1, col=2,
    )
    fig5.update_xaxes(title_text="Status", row=1, col=2)
    fig5.update_yaxes(title_text="Cell Count", row=1, col=2)

fig5.update_layout(title_text="Quality Metrics", height=400, showlegend=False)
figures.append(("Quality Metrics", fig5))

# ── Write HTML ─────────────────────────────────────────────────────────────────
import plotly.io as pio

html_parts = [
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>SPaLaTra Report</title>",
    "<style>body{font-family:sans-serif;margin:2em;} h1{color:#333;} h2{color:#555;border-bottom:1px solid #ccc;padding-bottom:4px;}</style>",
    "</head><body>",
    "<h1>SPaLaTra Annotation Report</h1>",
]

first = True
for title, fig in figures:
    html_parts.append(f"<h2>{title}</h2>")
    div_html = pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs="cdn" if first else False,
    )
    html_parts.append(div_html)
    first = False

html_parts.append("</body></html>")

os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
with open(args.output, "w", encoding="utf-8") as f:
    f.write("\n".join(html_parts))

print(f"Report written to {args.output}")
