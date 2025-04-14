#!/bin/bash

# Path to the Snakemake workflow file
WORKFLOW_FILE="brava_hits_common_condition_check.smk"

mkdir -p saige_outputs run_files

CORES=$(nproc)

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES --max-status-checks-per-second 0.01 \
    --keep-going --retries 3 --rerun-incomplete \
    -n \
    > "$LOGFILE" 2>&1

echo "Run complete. Log saved to $LOGFILE"
