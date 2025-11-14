#!/bin/bash

# Input and output variables
VCF="${1}" 
VARIANTS_COMMA="${2}"
MODELFILE="${3}"
VARIANCERATIO="${4}"
SPARSEGRM="${5}"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
CHR="${6}"
USE_NULL_VAR_RATIO="${7}"
P_T="${8}"

echo "null var vatio"
echo "${USE_NULL_VAR_RATIO}"

# Edge case
if [ $(bcftools index -n $VCF) -gt 0 ]; then 

  TMPFILE=$(mktemp)
  cmd=(step2_SPAtests.R
      --vcfFile=${VCF}
      --vcfFileIndex=${VCF}.csi
      --vcfField=GT
      --chrom="${CHR}"
      --minMAF=0
      --minMAC=10
      --GMMATmodelFile=${MODELFILE}
      --varianceRatioFile=${VARIANCERATIO}
      --LOCO=FALSE
      --is_Firth_beta=TRUE
      --pCutoffforFirth=0.10
      --is_output_moreDetails=TRUE
      --is_fastTest=TRUE
      --SAIGEOutputFile=${TMPFILE})

  if [ "${USE_NULL_VAR_RATIO,,}" = "false" ]; then
      cmd+=(--sparseGRMFile="${SPARSEGRM}"
      --sparseGRMSampleIDFile="${SPARSEGRMID}")
  fi

  # Run the command
  "${cmd[@]}"

  # Obtain the top hit as the first conditioning marker
  cond_M=$(sort -g -k13,13 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
  P_top=$(sort -g -k13,13 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $13}')

  echo "Lowest Pvalue in the sumstats file"
  echo "${cond_M}: ${P_top}"
  echo "Pvalue for significance for conditioning"
  echo ${P_T}

  # Add the top variant to the list of conditioning markers
  CONDITION=${cond_M}
  intFlag=$(python3 -c "print(1 if ${P_top} < ${P_T} else 0)")
  rm -f "${TMPFILE}"

  while [ "${intFlag}" -eq 1 ]
  do
    # Run the step2_SPAtests.R and redirect output to TMPFILE
    cond_cmd=("${cmd[@]}" --condition="${CONDITION}")
    "${cond_cmd[@]}"

    ncol=$(awk '{print NF; exit}' "${TMPFILE}")
    if [ "$ncol" -eq 29 ]; then
      cond_M=$(sort -g -k20,20 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
      P_top=$(sort -g -k20,20 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $20}')
    elif [ "$ncol" -eq 19 ]; then
      cond_M=$(sort -g -k18,18 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
      P_top=$(sort -g -k18,18 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $18}')
    else
      echo "Unexpected number of columns ($ncol) in $TMPFILE" >&2
      exit 1
    fi

    echo "Lowest Pvalue in the sumstats file"
    echo "${cond_M}: ${P_top}"
    echo "Pvalue for significance for conditioning"
    echo ${P_T}

    CONDITION_unordered="${CONDITION},${cond_M}"
    echo "conditioning..."
    
    # Write a small R script to ensure that the conditioning SNPs are 
    # in order
    CONDITION=$(Rscript scripts/sort_conditioning_snps.R --condition "${CONDITION_unordered}")
    echo $CONDITION

    intFlag=$(python3 -c "print(1 if ${P_top} < ${P_T} else 0)")
    rm -f "${TMPFILE}"
  done
else
  echo "No common variants present in the region" 
fi

# Note that the very last variant in "$CONDITION_unordered" is not significant any more after the final conditioning:
CONDITION=$(echo "$CONDITION_unordered" | awk -F',' 'NF>1 { for (i=1; i<NF; i++) printf "%s%s", $i, (i<NF-1 ? "," : "") }')
echo "${CONDITION}" > ${VARIANTS_COMMA}
