#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate brava_hits_common_condition_check

# Path to the Snakemake workflow file
cd ${SCRIPT_DIR}
# Make files present in same filepath, so that file locations
# in the list files are relative to the script directory
mv ${DATA_DIR} ${SCRIPT_DIR}
mv configs/${ANC}_config.yaml config.yaml
ls . 
echo "mv configs/${ANC}_config.yaml config.yaml"
WORKFLOW_FILE="brava_hits_common_stepwise_condition_check.smk"

mkdir -p saige_outputs run_files

CORES=$(nproc)

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES --max-status-checks-per-second 0.01 \
    --keep-going --retries 3 --rerun-incomplete --printshellcmds --verbose --forcerun all \
    > "$LOGFILE" 2>&1

echo "Run complete. Log saved to $LOGFILE"
cp ${LOGFILE} ${OUTPUT}/
