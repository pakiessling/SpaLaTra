# SPaLaTra — SPatial Label Transfer

SPaLaTra is a Snakemake pipeline that transfers cell type labels from a single-cell RNA-seq reference to spatial transcriptomics query datasets. It runs five independent annotation methods in parallel and combines them into a single plurality-consensus prediction, along with an interactive HTML report.

## Methods

| Method | Language | Key |
|---|---|---|
| [TACCO](https://github.com/simonwm/tacco) | Python | `tacco` |
| [SingleR](https://bioconductor.org/packages/SingleR/) | R | `singler` |
| [RCTD / spacexr](https://github.com/dmcable/spacexr) | R | `rctd` |
| [InSituType](https://github.com/Nanostring-Biostats/InSituType) | R | `insitutype` |
| [PhiSpace](https://github.com/jiadongm/PhiSpace) | R | `phispace` |

---

## Quick Start

```bash
# 1. Create and activate the conda environment
conda env create -f environment.yml
conda activate <env_name>

# 2. Edit the config file
#    Set input, output, ref, and optionally methods
nano config/config.yaml

# 3. Launch (SLURM cluster)
./launch.sh

# Dry-run to preview jobs without submitting
./launch.sh --dry-run

# Run locally (for testing, no SLURM)
./launch.sh --local
```

---

## Configuration

All settings are in **`config/config.yaml`**. Edit this file before running the pipeline.

```yaml
input:  "input/xenium_hs"          # directory of query .h5ad files
output: "result/xenium_hs"         # output directory
ref:    "references/scrna_hs.h5ad" # reference single-cell .h5ad

methods:          # remove entries to disable specific methods (min. 2)
  - tacco
  - singler
  - rctd
  - phispace
  - insitutype

embedding: spatial  # obsm key used for coordinates in the HTML report
```

To use an alternate config file (e.g., for a different experiment):

```bash
./launch.sh --config-file path/to/other_config.yaml
```

---

## Input Requirements

**Query `.h5ad` files** (all files in `input/` are processed):

| Requirement | Detail |
|---|---|
| `.X` | Raw counts — **not** log-normalised |
| `.obs.index` | Must be unique **across all query files**, not just within each |
| `.obsm['spatial']` | NumPy array with x/y coordinates (or set `embedding` to another key) |
| Layer names | Must not contain the string `"counts"` |

**Reference `.h5ad`**: single-cell dataset with raw counts in `.X` and cell type labels in `.obs`.

---

## Pipeline Architecture

```
config/config.yaml
        │
        ▼
Query .h5ad files + Reference .h5ad
        │
        ├── run_tacco.py        → tacco/{sample}_tacco.csv
        ├── run_singler.R       → singler/{sample}_singler.csv
        ├── run_rctd.R          → rctd/{sample}_rctd.csv
        ├── run_insitutype.R    → insitutype/{sample}_insitutype.csv
        └── run_phi_space.R     → phispace/{sample}_phispace.csv
                │
          combine.py  →  consensus.csv
                │
          report.py   →  report.html
```

Each method runs independently per sample. `combine.py` joins all per-cell predictions and computes a plurality consensus, using secondary-choice votes to break ties.

---

## Output

| File | Description |
|---|---|
| `{method}/{sample}_{method}.csv` | Per-cell predictions for each method and sample |
| `consensus.csv` | All methods joined; one row per cell, with `consensus`, `agreement_score`, and `is_ambiguous` columns |
| `report.html` | Interactive QC report (see below) |

### `consensus.csv` columns

| Column | Description |
|---|---|
| `tacco`, `singler`, `rctd`, `phispace`, `insitutype` | Per-method primary label |
| `tacco_2nd`, `phispace_2nd`, `rctd_2nd` | Secondary label (used for tie-breaking) |
| `singler_class`, `rctd_class` | Method-specific quality flags |
| `consensus` | Plurality-consensus label (`"unknown"` if unresolvable) |
| `agreement_score` | Fraction of methods agreeing with the consensus (0–1) |
| `is_ambiguous` | `True` when consensus is `"unknown"` |

### `report.html`

Self-contained interactive HTML with four sections:

1. **Per-method scatter plots** — cells coloured by label for each method and the consensus
2. **Pairwise agreement heatmap** — % of cells where each pair of methods agrees
3. **Label frequency bar chart** — cell type counts per method
4. **Quality metrics** — agreement score distribution and ambiguous cell count

---

## Disabling Methods

Remove entries from the `methods` list in `config/config.yaml`. At least two methods are required.

```yaml
# Example: run only three methods
methods:
  - tacco
  - singler
  - phispace
```

Only the selected methods are submitted to SLURM and included in the consensus.

---

## SLURM Configuration

SLURM settings live in **`slurm/config.yaml`** (the Snakemake workflow profile).

Default resources: **150 GB RAM, 5 CPUs, 12 h walltime**.

Per-rule overrides:

| Rule | RAM | CPUs |
|---|---|---|
| `singler` | 120 GB | 20 |
| `phispace` | 150 GB | 4 |
| `tacco` | 50 GB | 4 |
| `rctd` | 150 GB | 5 |
| `insitutype` | 50 GB | 4 |
| `consensus` | 50 GB | 1 |

> RCTD writes a placeholder CSV on error so downstream rules are never blocked.
