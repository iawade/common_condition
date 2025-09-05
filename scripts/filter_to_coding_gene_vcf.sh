#!/bin/bash

# Input parameters
INPUT_VCF="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
MAF_COMMON="$4" 
THREADS="$5"

EXPANDED_BED="run_files/bed/expanded_regions_${ENSEMBL_ID}.bed"

# Output files
OUTPUT_VCF="run_files/${ENSEMBL_ID}_${BP_DISTANCE}_${MAF_COMMON}.vcf.bgz"

# Use bcftools to filter VCF by the expanded BED regions and MAF threshold
# Using && which is the same as max(MAC > 40, MAF > $MAF_COMMON) ; unless I'm losing the plot
bcftools view --threads "$THREADS" -R "$EXPANDED_BED" "$INPUT_VCF" |
  bcftools filter --threads "$THREADS" -i "MAC > 40 && MAF > $MAF_COMMON" |
  bcftools view -Oz -o "${OUTPUT_VCF}"

bcftools index --csi -f "${OUTPUT_VCF}"

echo "Created indexed bgzipped VCF files for the gene ${ENSEMBL_ID}"
