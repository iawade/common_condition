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

echo -e "VCF=$VCF\nOUT=$OUT\nMIN_MAC=$MIN_MAC\nMODELFILE=$MODELFILE\nVARIANCERATIO=$VARIANCERATIO\nSPARSEGRM=$SPARSEGRM\nSPARSEGRMID=$SPARSEGRMID\nGROUPFILE=$GROUPFILE\nANNOTATIONS=$ANNOTATIONS\nCONDITION=$CONDITION"

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

# Append TMPFILE to OUT (header and all)
[[ -s "${TMPFILE}" ]] && cat "${TMPFILE}" >> "${OUT}"

# Also copy over single assoc files to verify case/control numbers
[[ -s "${TMPFILE}.singleAssoc.txt" ]] && cat "${TMPFILE}.singleAssoc.txt" >> "${OUT}.singleAssoc.txt"

# Ensure output exists for Snakemake
touch "${OUT}"

# Clean up
rm -f "${TMPFILE}"
rm -f "${TMPFILE}.singleAssoc.txt"
