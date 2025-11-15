#!/bin/bash

# Input and output variables
PLINK="${1}" 
VARIANTS_COMMA="${2}"
PHENOFILE="${3}"
COVARFILE="${4}"
PREDFILE="${5}"
CHR="${6}"
P_T="${7}"
PHENOCOL="${8}"
COVARCOLLIST="${9}"
CATEGCOVARCOLLIST="${10}"
CASE_CONTROL="${11}"

if [[ "${CASE_CONTROL}" == "binary" ]]; then
  trait_flag="--bt --firth --approx --pThresh 0.1"
elif [[ "${CASE_CONTROL}" == "quantitative" ]]; then
  trait_flag="--qt --apply-rint"
else
  echo "No trait type flag recognised" >&2
  exit 1
fi

# Edge case
if [ $(wc -l < ${PLINK}.bim) -gt 2 ]; then
  TMPFILE=$(mktemp)
  cmd=(regenie
      --step 2
      --bed ${PLINK}
      --phenoFile ${PHENOFILE}
      --covarFile ${COVARFILE}
      ${trait_flag}
      --pred ${PREDFILE}
      --minMAC 10
      --bsize 400
      --covarColList "${COVARCOLLIST}"
      --catCovarList "${CATEGCOVARCOLLIST}"
      --out ${TMPFILE}) # Add the covar cols etc

  # Run the command
  "${cmd[@]}"

  # Obtain the top hit as the first conditioning marker
  topline=$(sort -g -k12,12r "${TMPFILE}_${PHENOCOL}.regenie" | head -n 2 | tail -1)
  cond_M=$(awk '{print "chr"$1":"$2":"$4":"$5}' <<< "$topline")
  P_top=$(awk '{print $12}' <<< "$topline")

  echo "Lowest Pvalue in the sumstats file"
  echo "${cond_M}: ${P_top}"
  echo "Pvalue for significance for conditioning"
  echo ${P_T}
  P_T_log=$(awk -v v="${P_T}" 'BEGIN {print -log(v)/log(10)}')

  # Add the top variant to the list of conditioning markers, write to a tmp file
  echo ${cond_M} > ${TMPFILE}.cond.txt
  intFlag=$(python3 -c "print(1 if ${P_top} > ${P_T_log} else 0)")
  cond_cmd=("${cmd[@]}" --condition-list "${TMPFILE}.cond.txt")

  while [ "${intFlag}" -eq 1 ]; do
    "${cond_cmd[@]}"
    topline=$(sort -g -k12,12r "${TMPFILE}_${PHENOCOL}.regenie" | head -n 2 | tail -1)
    cond_M=$(awk '{print "chr"$1":"$2":"$4":"$5}' <<< "$topline")
    P_top=$(awk '{print $12}' <<< "$topline")
    echo "Lowest Pvalue in the sumstats file"
    echo "${cond_M}: ${P_top}"
    echo "Pvalue for significance for conditioning"
    echo ${P_T}
    intFlag=$(python3 -c "print(1 if ${P_top} > ${P_T_log} else 0)")
    echo ${cond_M} >> ${TMPFILE}.cond.txt
  done
else
  echo "No common variants present in the region" 
fi

sed '$d' ${TMPFILE}.cond.txt > ${VARIANTS_COMMA}
