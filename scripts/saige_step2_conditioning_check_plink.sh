#!/bin/bash

# Input and output variables
PLINK="${1}" 
OUT="${2}"
MIN_MAC="$3"
MODELFILE="$4"
VARIANCERATIO="$5"
SPARSEGRM="$6"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
GROUPFILE="$7"
ANNOTATIONS="$8"
CONDITION=$(cat "${9}")
MAX_MAF=${10}
USE_NULL_VAR_RATIO=${11}

echo -e "PLINK=$PLINK
    OUT=$OUT
    MIN_MAC=$MIN_MAC
    MODELFILE=$MODELFILE
    VARIANCERATIO=$VARIANCERATIO
    SPARSEGRM=$SPARSEGRM
    SPARSEGRMID=$SPARSEGRMID
    GROUPFILE=$GROUPFILE
    ANNOTATIONS=$ANNOTATIONS
    CONDITION=$CONDITION"

# /tmp is used to guard against weird edge cases if plink files are split in the
# middle of genes
TMPFILE=$(mktemp)

cmd=(step2_SPAtests.R
    --bimFile=${PLINK}.bim
    --bedFile=${PLINK}.bed
    --famFile=${PLINK}.fam
    --minMAF=0
    --minMAC=${MIN_MAC}
    --GMMATmodelFile=${MODELFILE}
    --varianceRatioFile=${VARIANCERATIO}
    --LOCO=FALSE
    --is_Firth_beta=TRUE
    --pCutoffforFirth=0.10
    --is_output_moreDetails=TRUE
    --is_fastTest=TRUE
    --SAIGEOutputFile=${TMPFILE}
    --groupFile=${GROUPFILE}
    --annotation_in_groupTest=${ANNOTATIONS}
    --is_output_markerList_in_groupTest=TRUE
    --is_single_in_groupTest=TRUE
    --maxMAF_in_groupTest=${MAX_MAF})

if [ -n "$CONDITION" ]; then
    cmd+=(--condition="${CONDITION}")
    if [ "${USE_NULL_VAR_RATIO,,}" = "false" ]; then
        cmd+=(--sparseGRMFile="${SPARSEGRM}"
        --sparseGRMSampleIDFile="${SPARSEGRMID}")
    fi
else
    echo "No conditioning required - no change in P-value"
fi

# Run the command
"${cmd[@]}"

# Append TMPFILE to OUT (header and all)
mv "${TMPFILE}" ${OUT}

# Also copy over single assoc files to verify case/control numbers
mv "${TMPFILE}.singleAssoc.txt" "${OUT}.singleAssoc.txt"
