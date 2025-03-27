#!/bin/bash

# Input and output variables
VCF="${1}" 
OUT="${2}"
MIN_MAC="$3"
MODELFILE="$4"
VARIANCERATIO="$5"
SPARSEGRM="$6"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
GROUPFILE="$7"
ANNOTATIONS="$8"
CONDITION=$(cat "${9}") 
# no error thrown if condition string is empty

TMPFILE=$(mktemp)

# Run the step2_SPAtests.R and redirect output to TMPFILE
step2_SPAtests.R \
        --vcfFile=${VCF} \
        --vcfFileIndex="${VCF}.csi" \
        --vcfField="GT" \
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
        --SAIGEOutputFile=${TMPFILE} \
        --groupFile=$GROUPFILE \
        --annotation_in_groupTest=$ANNOTATIONS \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --condition="$CONDITION" \
        --maxMissing=1

# Check if output file does not exist or if temporary file is empty
if [[ ! -e "${OUT}" || ! -s "${TMPFILE}" ]]; then
    echo "No variants found for the gene on chromosome ${CHR}. Exiting gracefully."
    touch "${OUT}"  # Ensure the file exists to prevent Snakemake errors
    rm ${TMPFILE}
    exit 0
fi

# If output file exists and temporary file is empty, append an empty file to OUT
if [[ -e "${OUT}" && ! -s "${TMPFILE}" ]]; then
    touch "${OUT}"
    rm ${TMPFILE}
    exit 0
fi

# Extract the header from the TMPFILE and append the rest of the content
HEADER=$(head -n 1 ${TMPFILE})
tail -n +2 ${TMPFILE} | sort -V -u >> ${OUT}
echo "${HEADER}" | cat - ${OUT} > temp && mv temp ${OUT}

rm ${TMPFILE}