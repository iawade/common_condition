#!/bin/bash

# Input and output variables
PLINK="${1}" 
VARIANTS_COMMA="${2}"
PHENOFILE="${3}"
COVARFILE="${4}"
PREDFILE="${5}" # This should be defined in at the start of the smk file
CHR="${6}"
P_T="${7}"
# CASE_CONTROL="${8}"

# conda activate regenie_env
# TMPFILE=$(mktemp)
# regenie \
#   --step 2 \
#   --bed run_files/ENSG00000167601_500000_0.0001 \
#   --phenoFile snakemake_eur/phenotypes/ukb.standing_height.20250508.tsv \
#   --covarFile snakemake_eur/covariates/ukb_brava_default_covariates.20250508.tsv \
#   --qt --apply-rint \
#   --pred run_files/Height_pred.list \
#   --minMAC 10 \
#   --bsize 400 \
#   --out ${TMPFILE}

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

  # # Obtain the top hit as the first conditioning marker
  # topline=$(sort -g -k13,13 "${TMPFILE}" | head -n 2 | tail -1)
  # cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  # P_top=$(awk '{print $13}' <<< "$topline")
  # lowest_pos_id=$(sort -g -k2,2 "${TMPFILE}" | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')

  # echo "Lowest Pvalue in the sumstats file"
  # echo "${cond_M}: ${P_top}"
  # echo "Pvalue for significance for conditioning"
  # echo ${P_T}
  # vars=2

  # while [ ${cond_M} = ${lowest_pos_id} ]; do
  #   echo "Weird edge case - SAIGE cannot condition on the first variant in the bim."
  #   vars=$((vars+1))
  #   topline=$(sort -g -k13,13 "${TMPFILE}" | head -n ${vars} | tail -1)
  #   cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  #   P_top=$(awk '{print $13}' <<< "$topline")
  # done
  # vars=2

  # # Add the top variant to the list of conditioning markers
  # CONDITION=${cond_M}
  # # intFlag=$(awk -v P_top="${P_top}" -v P_T="${P_T}" 'BEGIN{print (P_top<P_T)?1:0}')
  # intFlag=$(python3 -c "print(1 if ${P_top} < ${P_T} else 0)")

  # while [ "${intFlag}" -eq 1 ]
  # do
  #   # Run the step2_SPAtests.R and redirect output to TMPFILE
  #   cond_cmd=("${cmd[@]}" --condition="${CONDITION}")
  #   "${cond_cmd[@]}"

  #   ncol=$(awk '{print NF; exit}' "${TMPFILE}")
  #   if [ "$ncol" -eq 29 ]; then

  #     topline=$(sort -g -k20,20 "${TMPFILE}" | head -n 2 | tail -1)
  #     cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  #     P_top=$(awk '{print $20}' <<< "$topline")
  #     # Note that we need to re-determine the first variant in the output, as it seems to 
  #     # vary based on whether we are performing conditioning or not
  #     lowest_pos_id=$(sort -g -k2,2 "${TMPFILE}" | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
  #     echo "Lowest position variant's ID: ${lowest_pos_id}"

  #     while [ ${cond_M} = ${lowest_pos_id} ]; do
  #       echo "Weird edge case - SAIGE cannot condition on the first variant in the bim."
  #       vars=$((vars+1))
  #       topline=$(sort -g -k20,20 "${TMPFILE}" | head -n ${vars} | tail -1)
  #       cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  #       P_top=$(awk '{print $20}' <<< "$topline")
  #     done
  #     vars=2

  #   elif [ "$ncol" -eq 19 ]; then

  #     topline=$(sort -g -k18,18 "${TMPFILE}" | head -n 2 | tail -1)
  #     cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  #     P_top=$(awk '{print $18}' <<< "$topline")
  #     # Note that we need to re-determine the first variant in the output, as it seems to 
  #     # vary based on whether we are performing conditioning or not
  #     lowest_pos_id=$(sort -g -k2,2 "${TMPFILE}" | head -n 2 | tail -1 | awk '{print $1":"$2":"$4":"$5}')
  #     echo "Lowest position variant's ID: ${lowest_pos_id}"

  #     while [ ${cond_M} = ${lowest_pos_id} ]; do
  #       echo "Weird edge case - SAIGE cannot condition on the first variant in the bim."
  #       vars=$((vars+1))
  #       topline=$(sort -g -k18,18 "${TMPFILE}" | head -n ${vars} | tail -1)
  #       cond_M=$(awk '{print $1":"$2":"$4":"$5}' <<< "$topline")
  #       P_top=$(awk '{print $18}' <<< "$topline")
  #     done
  #     vars=2

  #   else
  #     echo "Unexpected number of columns ($ncol) in $TMPFILE" >&2
  #     exit 1
  #   fi

  #   echo "Lowest Pvalue in the sumstats file"
  #   echo "${cond_M}: ${P_top}"
  #   echo "Pvalue for significance for conditioning"
  #   echo ${P_T}

  #   CONDITION_unordered="${CONDITION},${cond_M}"
  #   echo "conditioning..."
    
  #   # Write a small R script to ensure that the conditioning SNPs are 
  #   # in order
  #   CONDITION=$(Rscript scripts/sort_conditioning_snps.R --condition "${CONDITION_unordered}")
  #   echo $CONDITION

  #   # intFlag=$(awk -v P_top="${P_top}" -v P_T="${P_T}" 'BEGIN{print (P_top<P_T)?1:0}')
  #   intFlag=$(python3 -c "print(1 if ${P_top} < ${P_T} else 0)")
  #   rm -f "${TMPFILE}"
  # done
else
  echo "No common variants present in the region" 
fi

# Note that the very last variant in "$CONDITION_unordered" is not significant any more after the final conditioning:
# CONDITION=$(echo "$CONDITION_unordered" | awk -F',' 'NF>1 { for (i=1; i<NF; i++) printf "%s%s", $i, (i<NF-1 ? "," : "") }')
CONDITION=""
echo "${CONDITION}" > ${VARIANTS_COMMA}

# # Move and compress the output
# mv *regenie ${OUT}.regenie
# gzip *regenie

# #rm tmp-*

# regenie requires info about whether it's cts or case/control - I can extract that from
# the phenotype json file
