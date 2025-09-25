#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate brava_hits_common_condition_check

# Path to the Snakemake workflow file
cd ${SCRIPT_DIR}
# Make files present in same filepath, so that file locations
# in the list files are relative to the script directory
mkdir ${DATA_DIR}

mv ${MODEL_DIR} ${DATA_DIR}
mv ${MATRIX_DIR} ${DATA_DIR}
mv ${VARIANCE_DIR} ${DATA_DIR}
mv ${GROUP_DIR} ${DATA_DIR}
VCF_DIR="${DATA_DIR}/vcf/combined"
mkdir -p ${VCF_DIR}
mv ${VCF} ${VCF_DIR}
mv ${VCF_CSI} ${VCF_DIR}

# move any already completed run information into the VM
mv ${RUN_DIR} final_run_files
mv ${SAIGE_DIR} final_saige_outputs

mv all-of-us_configs/${ANC}_config.yaml config.yaml
Rscript scripts/create_chr_specific_filepaths.r --chr ${CHR}

WORKFLOW_FILE="step2_final_group_tests.smk"

mkdir -p final_saige_outputs final_run_files/bed

CORES=$(nproc)

# Generate a timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="snakemake_run_${TIMESTAMP}.log"

# Run Snakemake with the specified options
echo "Starting a run of Snakemake workflow..."
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --touch
snakemake --snakefile "$WORKFLOW_FILE" --cores $CORES --jobs $CORES \
    --max-status-checks-per-second 0.01 --keep-going --rerun-incomplete \
    --printshellcmds --verbose --rerun-triggers code input params software-env \
     2>&1 | tee "$LOGFILE"

echo "Run complete. Log saved to $LOGFILE"
mv ${LOGFILE} ${OUTPUT}/
mv logs ${OUTPUT}/${ANC}_${CHR}_logs
mv final_saige_outputs ${OUTPUT}/${ANC}_${CHR}_final_saige_outputs
mv final_run_files ${OUTPUT}/${ANC}_${CHR}_final_run_files
mv brava_final_conditional_analysis_results.txt ${OUTPUT}/brava_${ANC}_${CHR}_final_conditional_analysis_results.txt
mv brava_final_conditional_analysis_results.txt.singleAssoc.txt ${OUTPUT}/brava_${ANC}_${CHR}_final_conditional_analysis_results.txt.singleAssoc.txt
