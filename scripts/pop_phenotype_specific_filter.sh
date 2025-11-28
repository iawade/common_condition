#!/bin/bash

# Input parameters
MODELFILE="$1" # Prevents any potential issues with gene symbols
CONDITION="$2" # File containing conditioning variants
CONDITION_BED="$3"
VCF="$4" # Input vcf
FILE="$5"

TMPFILE=$(mktemp)
Rscript scripts/filter_to_samples.R -m ${MODELFILE} -o ${TMPFILE}.sample

# Define a bed file first
plink2 --vcf $VCF --extract range ${CONDITION_BED} \
  --keep ${TMPFILE}.sample --mac 10 --make-bed --out ${TMPFILE} || true

if [[ -f ${TMPFILE}.bim ]]; then

    # Sort the .bim file variant ID
    awk 'BEGIN{OFS="\t"} { $2 = "chr" $1 ":" $4 ":" $6 ":" $5; print }' \
        ${TMPFILE}.bim > ${TMPFILE}.bim.tmp
    mv ${TMPFILE}.bim.tmp ${TMPFILE}.bim

    # Extract the exact variants
    plink2 --bfile ${TMPFILE} --extract $CONDITION --mac 10 \
        --make-bed --out ${TMPFILE}.tmp || true

    nvar=$(wc -l < ${TMPFILE}.tmp.bim)

    if [[ $nvar -eq 1 ]]; then
        cut -f2 ${TMPFILE}.tmp.bim > ${FILE}.txt
    elif [[ "$nvar" -gt 1 ]]; then
        plink2 --bfile ${TMPFILE}.tmp \
          --indep-pairwise 50 5 0.9 \
          --out ${FILE} || true
        # Finally, create a comma separated string from this
        if [[ -f ${FILE}.prune.in ]]; then
            paste -sd, ${FILE}.prune.in > ${FILE}.txt
        fi
    else
        # No conditioning variants present
        cat ${TMPFILE}.bim
        touch ${FILE}.txt
    fi
else
    # No conditioning variants present in the vcf file
    echo "None of the variants are present in the .vcf file"
    touch ${FILE}.txt
    nvar=$(wc -l < ${FILE}.txt)
    echo $nvar
fi
