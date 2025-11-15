#!/bin/bash

# Input and output variables
PLINK="${1}" 
OUT="${2}"
MIN_MAC="${3}"
PHENOFILE="${4}"
COVARFILE="${5}"
PHENOCOL=${6}
COVARCOLLIST=${7}
CATEGCOVARCOLLIST=${8}
PREDFILE="${9}"
ANNOTATIONFILE="${10}"
SETLISTFILE="${11}"
MASKDEF="${12}"
CONDITION="${13}"
MAX_MAF=${14}
CASE_CONTROL="${15}"

if [[ "${CASE_CONTROL}" == "binary" ]]; then
  trait_flag="--bt --firth --approx --pThresh 0.1"
elif [[ "${CASE_CONTROL}" == "quantitative" ]]; then
  trait_flag="--qt --apply-rint"
else
  echo "No trait type flag recognised" >&2
  exit 1
fi

# /tmp is used to guard against weird edge cases if plink files are split in the
# middle of genes
TMPFILE=$(mktemp)
cmd=(regenie
  --step 2
  --bed ${PLINK}
  --phenoFile ${PHENOFILE}
  --covarFile ${COVARFILE}
  ${trait_flag}
  --pred ${PREDFILE}
  --covarColList "${COVARCOLLIST}"
  --catCovarList "${CATEGCOVARCOLLIST}"
  --anno-file ${ANNOTATIONFILE}
  --set-list ${SETLISTFILE}
  --mask-def ${MASKFILE}
  --aaf-bins ${MAX_MAF}
  --vc-tests "skat,skato,acato"
  --minMAC ${MIN_MAC}
  --bsize 400
  --out ${TMPFILE})

if [ -s "$CONDITION" ]; then
    echo "There are variants to condition on"
    cmd+=(--condition-list "${CONDITION}")
fi

# Run the command
"${cmd[@]}"

# ALSO NEED TO RUN THE VARIANTS AS WELL (USING THE PREVIOUS COMMAND)
# Create the final output file for the conditioning variants
paste -sd, ${CONDITION} > "${TMPFILE}.conditioning.txt" && mv "${TMPFILE}.conditioning.txt" "${CONDITION}"
