#!/bin/bash

# Input and output variables
PLINK="${1}" 
VARIANTS_COMMA="${2}"
PHENOFILE="${3}"
COVARFILE="${4}"
PREDFILE="${5}" # This should be defined in at the start of the smk file
CHR="${6}"
P_T="${7}"
PHENOCOL=${8}
# CASE_CONTROL="${8}"

# conda activate regenie_env
# TMPFILE=$(mktemp)
# cmd=(regenie \
#   --step 2 \
#   --bed run_files/ENSG00000167601_500000_0.0001 \
#   --phenoFile snakemake_eur/phenotypes/ukb.standing_height.20250508.tsv \
#   --covarFile snakemake_eur/covariates/ukb_brava_default_covariates.20250508.tsv \
#   --qt --apply-rint \
#   --pred run_files/Height_pred.list \
#   --minMAC 10 \
#   --bsize 400 \
#   --out ${TMPFILE})
# "${cmd[@]}"
# # categorical covariate cols - use Nik's flags in the step 1 log he ran
# # phenotype_cols - don't need this, we let it be defined by the available loco files

# # We need to split the job up into tiny little jobs
# # only allowing gene-phenotype pairs

# # Decide what I need to move to the node and move it, and then start 
# # and interactive job. Make sure that we use the new docker.
# # then I can just iterate and get to the end

# # SAMPLE IDs to keep - equivalent as the mtx IDs I think
# # This should be a list containing the first two columns of the fam file - this should
# # be used instead of the sparse mtx sample IDs

# # Remove output prefix from .loco file
# # Generate this file on the fly - won't take long so just do it for all of the loco files
# # Define the pred.list for ourselves - we'll always just use a single trait
# # since we need to run for (gene, phenotype) pairs
# PRED="regenie_step1_${anc}_${PHENOCOL}_pred.list"
# LOCO="regenie_step1_${anc}_${PHENOCOL}_1.loco"
# PRED_LOCAL="$HOME/tmp-predfile.txt"
# cat ${PRED} | sed 's/\/home\/dnanexus\/out\/out\///g' > ${PRED_LOCAL} # Just rename the location of the predfile
# PRED=${PRED_LOCAL}
# head $PRED

# Determine whehter cts or binary
# If case_control is True 
#trait_flag="--bt --firth --approx --pThresh 0.1"
# else
trait_flag="--qt --apply-rint"

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

sed '$d' ${TMPFILE}.cond.txt | paste -sd, - > ${VARIANTS_COMMA}
