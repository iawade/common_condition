#!/bin/bash

# Path to the Snakemake workflow file
WORKFLOW_FILE="brava_hits_common_condition_check.smk"

mkdir -p saige_outputs

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores 1 --jobs 1 --max-status-checks-per-second 0.01 \
    2>&1 | tee snakemake_run.log 

echo "Run complete."	
