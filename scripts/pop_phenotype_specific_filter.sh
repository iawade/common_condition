#!/bin/bash

# Input parameters
MODELFILE="$1" # Prevents any potential issues with gene symbols
CONDITION="$2" # File containing conditioning variants
PLINK="$3" # Input plink fileset
FILE="$4"

TMPFILE=$(mktemp)
Rscript scripts/filter_to_samples.R -m ${MODELFILE} -o ${TMPFILE}.sample
plink2 --bfile ${PLINK} --extract ${CONDITION} --make-bed \
	--keep ${TMPFILE}.sample --mac 10 --make-bed \
	--out ${TMPFILE} || true
        
if [[ -f ${TMPFILE}.bim ]]; then
    nvar=$(wc -l < ${TMPFILE}.bim)
    if [[ $nvar -eq 1 ]]; then
        cut -f2 ${TMPFILE}.bim > ${FILE}.txt
    else
        plink2 --bfile ${TMPFILE} \
          --indep-pairwise 50 5 0.9 \
          --out ${FILE} || true
        # Finally, create a comma separated string from this
        if [[ -f ${FILE}.prune.in ]]; then
            paste -sd, ${FILE}.prune.in > ${FILE}.txt
        fi
    fi
else
    # No variants to prune, create an empty output
    echo "None of the variants are present in the plink .bim/.bed/.fam fileset"
    touch ${FILE}.txt
    nvar=$(wc -l < ${FILE}.txt)
    echo $nvar
fi
