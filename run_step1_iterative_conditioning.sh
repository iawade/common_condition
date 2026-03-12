#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate brava_hits_common_condition_check

mkdir -p saige_outputs run_files/bed

CORES=$(nproc)
full_mode=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--full)
      full_mode=true
      shift
      ;;
    *)
      shift
      ;;
  esac
done

if $full_mode; then
    echo "full GRM used to create sumstats"
    WORKFLOW_FILE="step1_iterative_conditioning_full_grm.smk"
else
    echo "sparse GRM used to create sumstats"
    WORKFLOW_FILE="step1_iterative_conditioning.smk"
fi

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES \
    --max-status-checks-per-second 0.01 --keep-going --rerun-incomplete \
    --printshellcmds --verbose --forcerun all \
    > "$LOGFILE" 2>&1

echo "Run complete. Log saved to $LOGFILE"
