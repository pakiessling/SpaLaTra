import argparse
import os

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
    df = df.idxmax(axis=1).to_frame()
    df.columns = ["phispace"]
    return df


def process_rctd(df):
    df = df[["spot_class", "first_type"]]
    df.columns = ["rctd_class", "rctd"]
    return df


def process_singler(df):
    df["singler_class"] = "good"
    df.loc[df["pruned.labels"].isna(), "singler_class"] = "bad"
    df = df.loc[:, ["labels", "singler_class"]]
    df.rename(columns={"labels": "singler"}, inplace=True)
    return df


def process_tacco(df):
    df = df.idxmax(axis=1).to_frame()
    df.columns = ["tacco"]
    return df


def process_insitutype(df):
    df.columns = ["insitutype"]
    return df


tacco = []
rctd = []
singler = []
phispace = []
insitutype = []

for root, dirs, files in os.walk(args.input):
    for file in files:
        if file.endswith(".csv"):
            path = os.path.join(root, file)
            df = pd.read_csv(path, index_col=0)

            if "phispace" in path:
                df = process_phispace(df)
                phispace.append(df)
            elif "rctd" in path:
                df = process_rctd(df)
                rctd.append(df)
            elif "singler" in path:
                df = process_singler(df)
                singler.append(df)
            elif "tacco" in path:
                df = process_tacco(df)
                tacco.append(df)
            elif "insitutype" in path:
                df = process_insitutype(df)
                insitutype.append(df)

            else:
                raise ValueError(f"Unknown file: {path}")

tacco = pd.concat(tacco)
rctd = pd.concat(rctd)
singler = pd.concat(singler)
phispace = pd.concat(phispace)
insitutype = pd.concat(insitutype)

combined = tacco.join(rctd).join(singler).join(phispace).join(insitutype)

combined["consensus"] = combined[["tacco", "rctd", "singler", "phispace", "insitutype"]].mode(
    axis=1
)[0]
combined.to_csv(args.output)
