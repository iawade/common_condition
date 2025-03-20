#!/bin/bash

# Set the maximum number of jobs
MAX_JOBS=$(nproc)

# Path to the Snakemake workflow file
WORKFLOW_FILE="brava_hits_common_condition_check.smk"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores "$MAX_JOBS" --jobs "$MAX_JOBS" --max-status-checks-per-second 0.01 \
    --keep-going \
    -n \
    2>&1 | tee snakemake_run.log 

echo "Run complete."	
