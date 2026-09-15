import argparse
import os
from collections import Counter

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "--input",
    type=str,
    help="Path to the input directory containing the CSV files (recursive).",
)
parser.add_argument(
    "--output",
    type=str,
    help="Path to the output CSV file where the combined results will be saved.",
)
parser.add_argument(
    "--methods",
    type=str,
    nargs="+",
    default=["tacco", "singler", "rctd", "phispace", "insitutype", "nnls", "tangram"],
    help="Active methods to include in the consensus.",
)
args = parser.parse_args()

ALL_METHODS = ["tacco", "singler", "rctd", "phispace", "insitutype", "nnls", "tangram"]
ACTIVE_METHODS = set(args.methods)


def process_argmax(df, name):
    """Top-1 and top-2 cell type per cell from a score matrix.

    Rows whose scores are all NaN (a method that failed on those cells) yield NaN
    rather than a bogus label; cells with a single finite score get no runner-up.
    """
    values = df.to_numpy(dtype=float)
    n_finite = np.isfinite(values).sum(axis=1)
    # -inf so NaN never wins, and so genuinely negative scores still compete
    # (PhiSpace scores are centred and can be negative).
    filled = np.where(np.isfinite(values), values, -np.inf)

    first = filled.argmax(axis=1)
    filled[np.arange(filled.shape[0]), first] = -np.inf
    second = filled.argmax(axis=1)

    cols = df.columns.to_numpy()
    primary = pd.Series(cols[first], index=df.index, dtype=object, name=name)
    secondary = pd.Series(
        cols[second], index=df.index, dtype=object, name=f"{name}_2nd"
    )
    primary[n_finite < 1] = np.nan
    secondary[n_finite < 2] = np.nan
    return pd.concat([primary, secondary], axis=1)


def process_phispace(df):
    return process_argmax(df, "phispace")


def process_rctd(df):
    result = df[["spot_class", "first_type"]].copy()
    result.columns = ["rctd_class", "rctd"]
    if "second_type" in df.columns:
        result["rctd_2nd"] = df["second_type"]
    else:
        result["rctd_2nd"] = "?"
    return result


def process_singler(df):
    out = df[["labels", "pruned.labels"]].copy()
    out["singler_class"] = "good"
    out.loc[out["pruned.labels"].isna(), "singler_class"] = "bad"
    out = out[["labels", "singler_class"]]
    out.rename(columns={"labels": "singler"}, inplace=True)
    return out


def process_tacco(df):
    return process_argmax(df, "tacco")


def process_nnls(df):
    return process_argmax(df, "nnls")


def process_tangram(df):
    return process_argmax(df, "tangram")


def process_insitutype(df):
    df = df.copy()
    df.columns = ["insitutype"]
    return df


PROCESSORS = {
    "phispace": process_phispace,
    "rctd": process_rctd,
    "singler": process_singler,
    "tacco": process_tacco,
    "insitutype": process_insitutype,
    "nnls": process_nnls,
    "tangram": process_tangram,
}


def plurality_consensus(row, primary_cols, secondary_cols):
    votes = [row[c] for c in primary_cols if pd.notna(row[c]) and row[c] != "?"]
    if not votes:
        return "unknown"
    counts = Counter(votes)
    max_count = max(counts.values())
    winners = [k for k, v in counts.items() if v == max_count]
    if len(winners) == 1:
        return winners[0]
    # Tie: add secondary votes to break it
    for col in secondary_cols:
        if col in row.index and pd.notna(row[col]) and row[col] != "?":
            votes.append(row[col])
    counts = Counter(votes)
    max_count = max(counts.values())
    winners = [k for k, v in counts.items() if v == max_count]
    return winners[0] if len(winners) == 1 else "unknown"


def classify(path):
    """Return (method, sample) for a result CSV, or (None, None) if unrecognised.

    Files are written as <output>/<method>/<sample>_<method>.csv, so the parent
    directory names the method and the basename names the sample.
    """
    method = os.path.basename(os.path.dirname(path))
    if method not in PROCESSORS:
        return None, None
    stem = os.path.basename(path)[: -len(".csv")]
    suffix = f"_{method}"
    sample = stem[: -len(suffix)] if stem.endswith(suffix) else stem
    return method, sample


# method -> list of (sample, path, frame); cells are keyed by (sample, cell) so that
# query datasets with colliding barcodes (e.g. Xenium ids, which restart at
# 'aaaaaaaa-1' in every run) stay distinct.
method_dfs = {method: [] for method in PROCESSORS}

for root, dirs, files in os.walk(args.input):
    dirs[:] = [d for d in dirs if d != "sections"]  # skip RCTD scatter intermediates
    for file in sorted(files):
        if not file.endswith(".csv"):
            continue
        path = os.path.join(root, file)
        method, sample = classify(path)
        if method is None:
            # e.g. a consensus.csv from a previous run sitting in the output root
            print(f"Skipping {path}: not under a <method>/ directory")
            continue
        if method not in ACTIVE_METHODS:
            continue
        df = pd.read_csv(path, index_col=0)
        df = PROCESSORS[method](df)
        df.index = pd.MultiIndex.from_arrays(
            [np.repeat(sample, len(df)), df.index.to_numpy()], names=["sample", "cell"]
        )
        method_dfs[method].append((sample, path, df))

for method in ACTIVE_METHODS:
    if not method_dfs[method]:
        raise ValueError(f"No {method} CSV files found in input directory")


# Validate cell uniqueness per method before concatenating. Duplicates are now only
# a problem within one sample -- across samples they are expected and handled.
def _check_unique(entries, method_name):
    seen = {}
    for sample, path, df in entries:
        dupes = df.index[df.index.duplicated()].get_level_values("cell").unique()
        if len(dupes) > 0:
            raise ValueError(f"Duplicate cell ids within {path}: {list(dupes[:10])}")
        if sample in seen:
            raise ValueError(
                f"Two {method_name} CSV files map to sample '{sample}': "
                f"{seen[sample]} and {path}"
            )
        seen[sample] = path


for method in ACTIVE_METHODS:
    _check_unique(method_dfs[method], method)

# Concatenate each method into a single frame
concatenated = {
    method: pd.concat([df for _, _, df in method_dfs[method]])
    for method in ACTIVE_METHODS
}

# Report index alignment across methods. Methods are not required to annotate the
# exact same set of cells: a method may legitimately drop cells (e.g. too few counts,
# per-section failures). Such cells become NaN on the outer join below and are simply
# ignored when computing the consensus. We build the union of all cells and warn about
# per-method coverage so genuine problems (wholesale drops, sample-name mismatches)
# are still visible in the logs.
full_index = pd.MultiIndex.from_arrays([[], []], names=["sample", "cell"])
for method in args.methods:
    full_index = full_index.union(concatenated[method].index)

for method in args.methods:
    missing = full_index.difference(concatenated[method].index)
    if len(missing) > 0:
        print(
            f"WARNING: '{method}' is missing {len(missing)} of {len(full_index)} cells "
            f"(e.g. {list(missing[:5])}); these will be NaN and ignored for that method."
        )

# Outer-join active methods so cells annotated by only some methods are retained.
combined = concatenated[args.methods[0]].reindex(full_index)
for method in args.methods[1:]:
    combined = combined.join(concatenated[method], how="outer")

# Build primary/secondary col lists from active methods
PRIMARY_COLS = [c for c in ALL_METHODS if c in ACTIVE_METHODS]
SECONDARY_COLS = [
    f"{c}_2nd"
    for c in ["tacco", "phispace", "rctd", "nnls", "tangram"]
    if c in ACTIVE_METHODS
]

# Normalize cell type names (some methods sanitize '/' → '_'; align all)
for col in PRIMARY_COLS + SECONDARY_COLS:
    if col in combined.columns:
        combined[col] = combined[col].str.replace("/", "_", regex=False)

combined["consensus"] = combined.apply(
    plurality_consensus,
    axis=1,
    primary_cols=PRIMARY_COLS,
    secondary_cols=SECONDARY_COLS,
)


# Quality metrics
def _agreement_score(row):
    valid = [row[c] for c in PRIMARY_COLS if pd.notna(row[c]) and row[c] != "?"]
    if not valid:
        return float("nan")
    return sum(v == row["consensus"] for v in valid) / len(valid)


combined["agreement_score"] = combined.apply(_agreement_score, axis=1)
combined["is_ambiguous"] = combined["consensus"] == "unknown"

# Write the original cell id as the index (as before) plus an explicit 'sample'
# column, so cells are still addressable when barcodes repeat across samples.
combined.insert(0, "sample", combined.index.get_level_values("sample"))
combined.index = combined.index.get_level_values("cell")
combined.index.name = "cell"
combined.to_csv(args.output)
