#!/bin/bash

# Path to the Snakemake workflow file
WORKFLOW_FILE="brava_hits_common_stepwise_condition_check.smk"

mkdir -p saige_outputs run_files

CORES=$(nproc)

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES --max-status-checks-per-second 0.01 \
    --keep-going --retries 3 --use-conda --conda-prefix "${HOME}/conda-envs" \
    --rerun-incomplete --printshellcmds --verbose --forcerun all \
    > "$LOGFILE" 2>&1

echo "Run complete. Log saved to $LOGFILE"
