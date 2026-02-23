# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**SPaLaTra** (SPatial Label Transfer) is a Snakemake pipeline that transfers cell type labels from single-cell RNA-seq reference datasets to spatial transcriptomics query datasets. It runs five independent annotation methods and combines them into a consensus prediction.

## Running the Pipeline

```bash
# Set up environment
conda env create -f enviroment.yml
conda activate <env_name>

# 1. Edit the config file
#    config/config.yaml  ← set input, output, ref, methods, etc.

# 2. Launch on SLURM
./launch.sh

# Or use an alternate config file
./launch.sh --config-file path/to/other_config.yaml

# Dry-run (preview jobs without submitting)
./launch.sh --dry-run

# Run locally (no SLURM)
./launch.sh --local
```

## Config File (`config/config.yaml`)

All pipeline settings live here. Edit this file before running.

| Key | Required | Description |
|---|---|---|
| `input` | yes | Directory of query `.h5ad` files (all processed) |
| `output` | yes | Output directory |
| `ref` | yes | Reference single-cell `.h5ad` file |
| `methods` | no | List of methods to run; defaults to all five. Minimum 2 required. |
| `embedding` | no | `obsm` key for report plots (default: `spatial`) |

Example — run only three methods:
```yaml
methods:
  - tacco
  - singler
  - phispace
```

## Architecture

The `Snakefile` orchestrates six rules: one per annotation method plus a final consensus step.

```
Query .h5ad files + Reference .h5ad
        │
        ├── run_tacco.py        (Python, TACCO)
        ├── run_singler.R       (R, SingleR)
        ├── run_rctd.R          (R, RCTD / spacexr)
        ├── run_insitutype.R    (R, InSituType)
        └── run_phi_space.R     (R, PhiSpace)
                │
          combine.py  →  consensus.csv
```

Each method script is invoked independently per sample and writes a CSV of per-cell predictions. `scripts/combine.py` joins all five CSVs and computes the mode across methods as the consensus label.

## Input Data Requirements

- `.X`: Raw counts (not log-normalized)
- `.obs.index`: Must be unique across **all** query datasets (not just within each file)
- `.obsm['spatial']`: NumPy array with x/y coordinates
- Layer names must not contain the string `"counts"`

## Per-Method Notes

| Method | Script | Output column used |
|---|---|---|
| TACCO | `run_tacco.py` | argmax of probability matrix |
| SingleR | `run_singler.R` | `labels`; quality flagged via `pruned.labels` |
| RCTD | `run_rctd.R` | `first_type` from `spot_class` |
| InSituType | `run_insitutype.R` | raw prediction |
| PhiSpace | `run_phi_space.R` | argmax of probability matrix |

RCTD writes a placeholder CSV on error so downstream rules are not blocked.

## SLURM Resource Allocation (`slurm/config.yaml`)

Default: 150 GB RAM, 5 CPUs, 12 h walltime. Per-rule overrides:

| Rule | RAM | CPUs |
|---|---|---|
| consensus | 50 GB | 1 |
| tacco | 50 GB | 4 |
| singler | 120 GB | 20 |
| insitutype | 50 GB | 4 |
| phispace | 150 GB | 4 |
