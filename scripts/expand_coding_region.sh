#!/bin/bash

# Input parameters
ENSEMBL_ID="$1" # Prevents any potential issues with gene symbols
BP_DISTANCE="$2" # TODO Error handling / kb vs just the number

# Output files
EXPANDED_BED="run_files/bed/expanded_regions_${ENSEMBL_ID}.bed"
EXPANDED_CODING_BED="run_files/bed/expanded_coding_regions_${ENSEMBL_ID}.bed"

awk -v BP_DISTANCE="$BP_DISTANCE" 'BEGIN {OFS="\t"} {
    start = $2 - BP_DISTANCE;
    end = $3 + BP_DISTANCE;

    if (start < 0) start = 0;
    
    # Add "chr" prefix if chromosome is numeric
    if ($1 ~ /^[0-9XY]+$/) $1 = "chr" $1;

    # Print the entire expanded region
    print $1, start, end;
}' "run_files/bed/${ENSEMBL_ID}.bed" > ${EXPANDED_BED}

bedtools intersect -a ${EXPANDED_BED} -b data/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed > "${EXPANDED_CODING_BED}"

echo "Created ${EXPANDED_BED} and ${EXPANDED_CODING_BED}, a padded coding region around the gene"
