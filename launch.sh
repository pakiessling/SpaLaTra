#!/usr/bin/env bash
# launch.sh — Start the SPaLaTra pipeline on a SLURM HPC cluster.
#
# Pipeline settings are read from config/config.yaml by default.
# Edit that file to set input, output, ref, methods, etc.
#
# Usage:
#   ./launch.sh                                    # use config/config.yaml
#   ./launch.sh --config-file path/to/config.yaml  # use alternate config
#   ./launch.sh --dry-run                          # preview jobs without submitting
#   ./launch.sh --local                            # run locally (no SLURM)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
CONFIG_FILE="$SCRIPT_DIR/config/config.yaml"
DRY_RUN=""
LOCAL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config-file) CONFIG_FILE="$2"; shift 2 ;;
        --dry-run)     DRY_RUN="--dry-run"; shift ;;
        --local)       LOCAL="1"; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

if ! command -v snakemake &>/dev/null; then
    echo "ERROR: snakemake not found in PATH." >&2
    echo "Activate the conda environment that has Snakemake installed." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build snakemake command
# ---------------------------------------------------------------------------
if [[ -n "$LOCAL" ]]; then
    CMD="snakemake \
        --snakefile \"$SCRIPT_DIR/Snakefile\" \
        --configfile \"$CONFIG_FILE\" \
        --conda-frontend conda \
        --use-conda \
        --cores 4 \
        $DRY_RUN"
else
    CMD="snakemake \
        --snakefile \"$SCRIPT_DIR/Snakefile\" \
        --configfile \"$CONFIG_FILE\" \
        --executor slurm \
        --conda-frontend conda \
        --workflow-profile \"$SCRIPT_DIR/slurm\" \
        --rerun-incomplete \
        $DRY_RUN"
fi

echo "Config: $CONFIG_FILE"
echo "Running: $CMD"
echo ""

if [[ -n "$DRY_RUN" ]]; then
    eval "$CMD"
else
    mkdir -p "$SCRIPT_DIR/logs"
    LOG="$SCRIPT_DIR/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
    nohup bash -c "$CMD" > "$LOG" 2>&1 &
    echo "Pipeline started in background. PID: $!"
    echo "Monitor progress: tail -f $LOG"
fi
