#!/bin/bash

# Input and output variables
VCF="${1}" 
VARIANTS_COMMA="${2}"
MODELFILE="${3}"
VARIANCERATIO="${4}"
SPARSEGRM="${5}"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
CHR="${6}"

TMPFILE=$(mktemp)

# Run the step2_SPAtests.R and redirect output to TMPFILE
step2_SPAtests.R \
        --vcfFile="${VCF}" \
        --vcfFileIndex="${VCF}.csi" \
        --vcfField="GT" \
        --chrom="${CHR}" \
        --minMAF=0 \
        --minMAC=10 \
        --GMMATmodelFile="${MODELFILE}" \
        --varianceRatioFile="${VARIANCERATIO}" \
        --sparseGRMFile="${SPARSEGRM}" \
        --sparseGRMSampleIDFile="${SPARSEGRMID}" \
        --LOCO=FALSE \
        --pCutoffforFirth=0.01 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --SAIGEOutputFile="${TMPFILE}"

# Signficant threshold for p-values
P_T=5e-3

# Obtain the top hit as the first conditioning marker
cond_M=$(sort -g -k13,13 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
P_top=$(sort -g -k13,13 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $13}')

echo 'Lowest Pvalue in the sumstats file'
echo ${P_top}
echo "Pvalue for significance for conditioning"
echo ${P_T}

# Add the top variant to the list of conditioning markers
CONDITION=${cond_M}
intFlag=$(awk -v P_top="${P_top}" -v P_T="${P_T}" 'BEGIN{print (P_top<P_T)?1:0}')
rm -f "${TMPFILE}"

while [ "${intFlag}" -eq 1 ]
do
  # Run the step2_SPAtests.R and redirect output to TMPFILE
  # TMPFILE=$(mktemp)

  step2_SPAtests.R \
          --vcfFile="${VCF}" \
          --vcfFileIndex="${VCF}.csi" \
          --vcfField="GT" \
          --chrom="${CHR}" \
          --minMAF=0 \
          --minMAC=10 \
          --GMMATmodelFile="${MODELFILE}" \
          --varianceRatioFile="${VARIANCERATIO}" \
          --sparseGRMFile="${SPARSEGRM}" \
          --sparseGRMSampleIDFile="${SPARSEGRMID}" \
          --LOCO=FALSE \
          --pCutoffforFirth=0.01 \
          --is_output_moreDetails=TRUE \
          --is_fastTest=TRUE \
          --SAIGEOutputFile="${TMPFILE}" \
          --condition="${CONDITION}"

  ncol=$(awk '{print NF; exit}' "${TMPFILE}")

  if [ "$ncol" -eq 29 ]; then
    cond_M=$(sort -g -k20,20 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
    P_top=$(sort -g -k20,20 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $20}')
  elif [ "$ncol" -eq 19 ]; then
    cond_M=$(sort -g -k18,18 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
    P_top=$(sort -g -k18,18 ${TMPFILE} | head -n 2 | tail -1 | awk '{print $18}')
  else
    echo "Unexpected number of columns ($ncol) in $TMPFILE" #>&2
    # exit 1
  fi

  echo 'Lowest Pvalue in the sumstats file'
  echo ${P_top}
  echo "Pvalue for significance for conditioning"
  echo ${P_T}

  CONDITION="${CONDITION},${cond_M}"
  echo "conditioning..."
  # Write a small R script to ensure that the conditioning SNPs are 
  # in order
  CONDITION=$(Rscript scripts/sort_conditioning_snps.R --condition "${CONDITION}")
  echo $CONDITION

  intFlag=$(awk -v P_top="${P_top}" -v P_T="${P_T}" 'BEGIN{print (P_top<P_T)?1:0}')
  rm -f "${TMPFILE}"
done

# Note that the very last variant in "$CONDITION" is not significant any more after the final conditioning:
CONDITION=$(echo "$CONDITION" | awk -F',' 'NF>1 { for (i=1; i<NF; i++) printf "%s%s", $i, (i<NF-1 ? "," : "") }')
echo "${CONDITION}" > ${VARIANTS_COMMA}
