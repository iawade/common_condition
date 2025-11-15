#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate brava_hits_common_condition_check

# Make files present in same filepath, so that file locations
# in the list files are relative to the script directory
mv uk-biobank_configs/${ANC}_config_regenie.yaml config.yaml
Rscript scripts/create_chr_specific_filepaths.r --chr ${CHR} -p ${P_T}

# Path to the Snakemake workflow file
WORKFLOW_FILE="step1_iterative_conditioning_regenie.smk"
CORES=$(nproc)

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Tidy up (just in case)
rm brava_stepwise_conditional_analysis_results*.txt

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --touch
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES --latency-wait 60 \
    --max-status-checks-per-second 0.01 \
    --printshellcmds --verbose --rerun-triggers code input params software-env \
     2>&1 | tee "$LOGFILE"

echo "Run complete. Log saved to $LOGFILE"
