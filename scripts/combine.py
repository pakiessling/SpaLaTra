import argparse
import os
from collections import Counter

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
args = parser.parse_args()


def process_phispace(df):
    # Primary: argmax col
    primary = df.idxmax(axis=1).to_frame(name="phispace")
    # Secondary: argmax after zeroing the primary column
    df2 = df.copy()
    for idx in df2.index:
        df2.loc[idx, primary.loc[idx, "phispace"]] = 0
    secondary = df2.idxmax(axis=1).to_frame(name="phispace_2nd")
    return primary.join(secondary)


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
    # Primary: argmax col
    primary = df.idxmax(axis=1).to_frame(name="tacco")
    # Secondary: argmax after zeroing the primary column
    df2 = df.copy()
    for idx in df2.index:
        df2.loc[idx, primary.loc[idx, "tacco"]] = 0
    secondary = df2.idxmax(axis=1).to_frame(name="tacco_2nd")
    return primary.join(secondary)


def process_insitutype(df):
    df = df.copy()
    df.columns = ["insitutype"]
    return df


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


tacco_dfs = []
rctd_dfs = []
singler_dfs = []
phispace_dfs = []
insitutype_dfs = []

for root, dirs, files in os.walk(args.input):
    for file in sorted(files):
        if file.endswith(".csv"):
            path = os.path.join(root, file)
            df = pd.read_csv(path, index_col=0)

            if "phispace" in path:
                df = process_phispace(df)
                phispace_dfs.append(df)
            elif "rctd" in path:
                df = process_rctd(df)
                rctd_dfs.append(df)
            elif "singler" in path:
                df = process_singler(df)
                singler_dfs.append(df)
            elif "tacco" in path:
                df = process_tacco(df)
                tacco_dfs.append(df)
            elif "insitutype" in path:
                df = process_insitutype(df)
                insitutype_dfs.append(df)
            else:
                raise ValueError(f"Unknown file: {path}")

# Empty-list guards
if not tacco_dfs:
    raise ValueError("No tacco CSV files found in input directory")
if not rctd_dfs:
    raise ValueError("No rctd CSV files found in input directory")
if not singler_dfs:
    raise ValueError("No singler CSV files found in input directory")
if not phispace_dfs:
    raise ValueError("No phispace CSV files found in input directory")
if not insitutype_dfs:
    raise ValueError("No insitutype CSV files found in input directory")

# Validate index uniqueness per method before concatenating
def _check_unique(dfs, method_name):
    all_indices = pd.Index([idx for df in dfs for idx in df.index])
    dupes = all_indices[all_indices.duplicated()].unique()
    if len(dupes) > 0:
        raise ValueError(
            f"Duplicate cell indices found across {method_name} CSV files: {list(dupes[:10])}"
        )

_check_unique(tacco_dfs, "tacco")
_check_unique(rctd_dfs, "rctd")
_check_unique(singler_dfs, "singler")
_check_unique(phispace_dfs, "phispace")
_check_unique(insitutype_dfs, "insitutype")

tacco = pd.concat(tacco_dfs)
rctd = pd.concat(rctd_dfs)
singler = pd.concat(singler_dfs)
phispace = pd.concat(phispace_dfs)
insitutype = pd.concat(insitutype_dfs)

combined = tacco.join(rctd).join(singler).join(phispace).join(insitutype)

PRIMARY_COLS = ["tacco", "rctd", "singler", "phispace", "insitutype"]
SECONDARY_COLS = ["tacco_2nd", "phispace_2nd", "rctd_2nd"]

combined["consensus"] = combined.apply(
    plurality_consensus, axis=1, primary_cols=PRIMARY_COLS, secondary_cols=SECONDARY_COLS
)

# Quality metrics
def _agreement_score(row):
    valid = [row[c] for c in PRIMARY_COLS if pd.notna(row[c]) and row[c] != "?"]
    if not valid:
        return float("nan")
    return sum(v == row["consensus"] for v in valid) / len(valid)

combined["agreement_score"] = combined.apply(_agreement_score, axis=1)
combined["is_ambiguous"] = combined["consensus"] == "unknown"

combined.to_csv(args.output)
