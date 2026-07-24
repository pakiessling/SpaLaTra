import argparse
import glob
import os

import anndata
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.metrics import cohen_kappa_score, f1_score

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

# ── Section 1b: Pairwise Cohen's kappa & macro-F1 ─────────────────────────────
# Raw agreement is inflated by imbalanced label frequencies (two methods that both
# over-call a dominant type look concordant). Cohen's kappa corrects for chance
# agreement, and macro-F1 weights every cell type equally so a method that disagrees
# on rare types is penalised. A method that scores well on % agreement but poorly on
# kappa/macro-F1 (RCTD is a common culprit) is a candidate for exclusion.
def _pairwise_matrix(metric_fn):
    mat = np.full((n, n), np.nan)
    for i, m1 in enumerate(PRIMARY_METHODS):
        for j, m2 in enumerate(PRIMARY_METHODS):
            if i == j:
                mat[i, j] = 1.0  # a method vs itself: perfect κ / F1 by definition
                continue
            valid = consensus[[m1, m2]].dropna()
            valid = valid[(valid[m1] != "?") & (valid[m2] != "?")]
            if len(valid) == 0:
                continue
            mat[i, j] = metric_fn(valid[m1], valid[m2])
    return mat

if n >= 2:
    kappa_mat = _pairwise_matrix(lambda a, b: cohen_kappa_score(a, b))
    kappa_df = pd.DataFrame(kappa_mat, index=PRIMARY_METHODS, columns=PRIMARY_METHODS)
    fig_kappa = px.imshow(
        kappa_df,
        text_auto=".2f",
        color_continuous_scale="RdBu",
        zmin=-1,
        zmax=1,
        title="Pairwise Cohen's κ (chance-corrected agreement; 1 = perfect, 0 = chance)",
        labels=dict(color="Cohen's κ"),
    )
    fig_kappa.update_layout(height=500)
    figures.append(("Pairwise Cohen's κ", fig_kappa))

    f1_mat = _pairwise_matrix(
        lambda a, b: f1_score(a, b, average="macro", zero_division=0)
    )
    f1_df = pd.DataFrame(f1_mat, index=PRIMARY_METHODS, columns=PRIMARY_METHODS)
    fig_f1 = px.imshow(
        f1_df,
        text_auto=".2f",
        color_continuous_scale="Greens",
        zmin=0,
        zmax=1,
        title="Pairwise macro-F1 (row = reference labels, per-cell-type mean)",
        labels=dict(color="macro-F1"),
    )
    fig_f1.update_layout(height=500)
    figures.append(("Pairwise macro-F1", fig_f1))

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

# ── Section 4: Cell type confusion heatmap ────────────────────────────────────
if "consensus" in consensus.columns and PRIMARY_METHODS:
    # For each cell, record every (consensus_label, method_label) pair where they disagree
    valid_consensus = consensus[consensus["consensus"].notna() & (consensus["consensus"] != "unknown")]
    all_labels = sorted(
        set(valid_consensus["consensus"].unique())
        | set(valid_consensus[PRIMARY_METHODS].stack().unique())
    )
    label_index = {l: i for i, l in enumerate(all_labels)}
    n_labels = len(all_labels)
    confusion = np.zeros((n_labels, n_labels), dtype=int)

    for _, row in valid_consensus.iterrows():
        c = row["consensus"]
        if c not in label_index:
            continue
        ci = label_index[c]
        for m in PRIMARY_METHODS:
            ml = row[m]
            if pd.isna(ml) or ml == "?" or ml == c:
                continue
            if ml in label_index:
                confusion[ci, label_index[ml]] += 1

    # Row-normalize: % of consensus-X cells called Y by at least one method
    row_sums = valid_consensus["consensus"].value_counts()
    norm = np.zeros_like(confusion, dtype=float)
    for i, label in enumerate(all_labels):
        total = row_sums.get(label, 0)
        if total > 0:
            norm[i, :] = confusion[i, :] / total * 100

    confusion_df = pd.DataFrame(norm, index=all_labels, columns=all_labels)
    fig4 = px.imshow(
        confusion_df,
        text_auto=".0f",
        color_continuous_scale="Reds",
        title="Cell Type Confusion (% of consensus-X cells called Y by ≥1 method)",
        labels=dict(x="Called by method", y="Consensus label", color="% cells"),
        aspect="auto",
    )
    fig4.update_layout(height=max(400, 25 * n_labels))
    figures.append(("Cell Type Confusion", fig4))

# ── Section 5: Label frequency bar chart ──────────────────────────────────────
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

# ── Section 6: Quality metrics ────────────────────────────────────────────────
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
