#!/bin/bash

# Input and output variables
VCF="${1}" 
OUT="${2}"
CHR="$3"
MIN_MAC="$4"
MODELFILE="$5"
VARIANCERATIO="$6"
SPARSEGRM="$7"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
GROUPFILE="$8"
ANNOTATIONS="$9"
CONDITION=$(cat "${10}") 
# no error thrown if condition string is empty

step2_SPAtests.R \
        --vcfFile=${VCF} \
        --vcfFileIndex="${VCF}.csi"\
        --vcfField="GT" \
        --chrom="$CHR" \
        --minMAF=0 \
        --minMAC=${MIN_MAC} \
        --GMMATmodelFile=${MODELFILE} \
        --varianceRatioFile=${VARIANCERATIO} \
        --sparseGRMFile=${SPARSEGRM} \
        --sparseGRMSampleIDFile=${SPARSEGRMID} \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --SPAcutoff=0.5 \
        --pCutoffforFirth=0.10 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${OUT} \
        --groupFile=$GROUPFILE \
        --annotation_in_groupTest=$ANNOTATIONS \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --condition="$CONDITION" \
        --maxMissing=1 \
        #--impute_method="mean"

# Check if output file is empty
if [[ ! -s "${OUT}" ]]; then
    echo "No variants found for the gene on chromosome ${CHR}. Exiting gracefully."
    touch "${OUT}"  # Ensure the file exists to prevent Snakemake errors
    exit 0
fi