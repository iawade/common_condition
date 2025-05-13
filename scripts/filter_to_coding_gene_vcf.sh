#!/bin/bash

# Input parameters
INPUT_VCF="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
MAF_COMMON="$4" 
THREADS="$5"

# Output files
OUTPUT_VCF="run_files/${ENSEMBL_ID}_${BP_DISTANCE}_${MAF_COMMON}.vcf.bgz"

# Expand the BED regions for query and filter to coding regions
## include common variation within the gene of interest too 
EXPANDED_BED="run_files/expanded_regions_${ENSEMBL_ID}.bed"
awk -v BP_DISTANCE="$BP_DISTANCE" 'BEGIN {OFS="\t"} {
    start = $2 - BP_DISTANCE;
    end = $3 + BP_DISTANCE;

    if (start < 0) start = 0;
    
    # Add "chr" prefix if chromosome is numeric
    if ($1 ~ /^[0-9X]+$/) $1 = "chr" $1;

    # Print the entire expanded region
    print $1, start, end;
}' "run_files/${ENSEMBL_ID}.bed" | bedtools intersect -a stdin -b protein_coding_regions_hg38_no_padding_no_UTR_v39.bed > "$EXPANDED_BED"

# Use bcftools to filter VCF by the expanded BED regions and MAF threshold
# Using && which is the same as max(MAC > 40, MAF > $MAF_COMMON) ; unless I'm losing the plot
bcftools view --threads "$THREADS" -R "$EXPANDED_BED" "$INPUT_VCF" |
  bcftools filter --threads "$THREADS" -i "MAC > 40 && MAF > $MAF_COMMON" |
  bcftools view -Oz -o "${OUTPUT_VCF}"

bcftools index --csi -f "${OUTPUT_VCF}"

echo "Created indexed bgzipped VCF files for the gene ${ENSEMBL_ID}"
