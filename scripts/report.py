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
parser.add_argument("--input", required=True, help="Directory of query .h5ad files")
parser.add_argument("--output", required=True, help="Output HTML path")
parser.add_argument(
    "--embedding", default="spatial", help="obsm key for spatial coordinates"
)
args = parser.parse_args()

# ── Load spatial coordinates ──────────────────────────────────────────────────
h5ad_files = glob.glob(os.path.join(args.input, "*.h5ad"))
if not h5ad_files:
    raise FileNotFoundError(f"No .h5ad files found in {args.input}")

coords_parts = []
for path in h5ad_files:
    adata = anndata.read_h5ad(path)
    if args.embedding not in adata.obsm:
        raise KeyError(
            f"Embedding key '{args.embedding}' not found in {path}. "
            f"Available keys: {list(adata.obsm.keys())}"
        )
    xy = adata.obsm[args.embedding]
    coords_parts.append(
        pd.DataFrame(
            {"x": xy[:, 0], "y": xy[:, 1]},
            index=adata.obs_names,
        )
    )
coords_df = pd.concat(coords_parts)

# ── Load consensus ────────────────────────────────────────────────────────────
consensus = pd.read_csv(args.consensus, index_col=0)

# Join coordinates to consensus
data = coords_df.join(consensus, how="inner")

METHOD_COLS = ["tacco", "singler", "rctd", "phispace", "insitutype", "consensus"]
PRESENT_METHODS = [c for c in METHOD_COLS if c in data.columns]

# Unified color mapping across all labels
all_labels = pd.unique(data[PRESENT_METHODS].values.ravel())
all_labels = [l for l in all_labels if pd.notna(l)]
palette = px.colors.qualitative.Alphabet
color_map = {
    label: palette[i % len(palette)] for i, label in enumerate(sorted(all_labels))
}

figures = []

# ── Section 1: Per-method scatter plots ───────────────────────────────────────
n_methods = len(PRESENT_METHODS)
cols = 3
rows = (n_methods + cols - 1) // cols

fig1 = make_subplots(
    rows=rows,
    cols=cols,
    subplot_titles=PRESENT_METHODS,
    shared_xaxes=False,
    shared_yaxes=False,
)

seen_labels = set()
for i, method in enumerate(PRESENT_METHODS):
    row = i // cols + 1
    col = i % cols + 1
    sub = data[["x", "y", method]].dropna()
    for label in sub[method].unique():
        mask = sub[method] == label
        show_legend = label not in seen_labels
        seen_labels.add(label)
        fig1.add_trace(
            go.Scatter(
                x=sub.loc[mask, "x"],
                y=sub.loc[mask, "y"],
                mode="markers",
                marker=dict(size=3, color=color_map.get(label, "#888")),
                name=label,
                legendgroup=label,
                showlegend=show_legend,
                text=label,
                hovertemplate=f"{method}: {label}<extra></extra>",
            ),
            row=row,
            col=col,
        )
    fig1.update_xaxes(showticklabels=False, row=row, col=col)
    fig1.update_yaxes(showticklabels=False, scaleanchor=f"x{i+1}" if i > 0 else "x", row=row, col=col)

fig1.update_layout(
    title_text="Per-Method Cell Type Annotations",
    height=400 * rows,
    legend=dict(itemsizing="constant"),
)
figures.append(("Per-Method Scatter Plots", fig1))

# ── Section 2: Pairwise agreement heatmap ────────────────────────────────────
primary_methods = [c for c in ["tacco", "singler", "rctd", "phispace", "insitutype"] if c in data.columns]
n = len(primary_methods)
agreement_matrix = np.zeros((n, n))
for i, m1 in enumerate(primary_methods):
    for j, m2 in enumerate(primary_methods):
        valid = data[[m1, m2]].dropna()
        valid = valid[(valid[m1] != "?") & (valid[m2] != "?")]
        if len(valid) == 0:
            agreement_matrix[i, j] = float("nan")
        else:
            agreement_matrix[i, j] = (valid[m1] == valid[m2]).mean() * 100

agreement_df = pd.DataFrame(agreement_matrix, index=primary_methods, columns=primary_methods)
fig2 = px.imshow(
    agreement_df,
    text_auto=".1f",
    color_continuous_scale="Blues",
    zmin=0,
    zmax=100,
    title="Pairwise Method Agreement (% cells with same label)",
    labels=dict(color="Agreement (%)"),
)
fig2.update_layout(height=500)
figures.append(("Pairwise Agreement Heatmap", fig2))

# ── Section 3: Label frequency bar chart ──────────────────────────────────────
freq_rows = []
for method in PRESENT_METHODS:
    counts = data[method].value_counts()
    for label, count in counts.items():
        freq_rows.append({"Method": method, "Label": label, "Count": int(count)})
freq_df = pd.DataFrame(freq_rows)

fig3 = px.bar(
    freq_df,
    x="Label",
    y="Count",
    color="Method",
    barmode="group",
    title="Cell Type Label Frequency per Method",
)
fig3.update_layout(height=500, xaxis_tickangle=-45)
figures.append(("Label Frequency", fig3))

# ── Section 4: Quality metrics ────────────────────────────────────────────────
fig4 = make_subplots(rows=1, cols=2, subplot_titles=["Agreement Score Distribution", "Ambiguous Cell Count"])

if "agreement_score" in data.columns:
    scores = data["agreement_score"].dropna()
    fig4.add_trace(
        go.Histogram(x=scores, nbinsx=20, name="Agreement Score", marker_color="#636EFA"),
        row=1, col=1,
    )
    fig4.update_xaxes(title_text="Agreement Score (0–1)", row=1, col=1)
    fig4.update_yaxes(title_text="Cell Count", row=1, col=1)

if "is_ambiguous" in data.columns:
    ambig_counts = data["is_ambiguous"].value_counts().reset_index()
    ambig_counts.columns = ["Ambiguous", "Count"]
    ambig_counts["Ambiguous"] = ambig_counts["Ambiguous"].map({True: "Ambiguous", False: "Resolved"})
    fig4.add_trace(
        go.Bar(
            x=ambig_counts["Ambiguous"],
            y=ambig_counts["Count"],
            name="Ambiguity",
            marker_color=["#EF553B", "#00CC96"],
        ),
        row=1, col=2,
    )
    fig4.update_xaxes(title_text="Status", row=1, col=2)
    fig4.update_yaxes(title_text="Cell Count", row=1, col=2)

fig4.update_layout(title_text="Quality Metrics", height=400, showlegend=False)
figures.append(("Quality Metrics", fig4))

# ── Write HTML ────────────────────────────────────────────────────────────────
html_parts = [
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>SPaLaTra Report</title>",
    "<style>body{font-family:sans-serif;margin:2em;} h1{color:#333;} h2{color:#555;border-bottom:1px solid #ccc;padding-bottom:4px;}</style>",
    "</head><body>",
    "<h1>SPaLaTra Annotation Report</h1>",
]

import plotly.io as pio

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
