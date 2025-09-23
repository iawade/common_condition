#!/bin/bash

# Input parameters
INPUT_VCF="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
THREADS="$4"

EXPANDED_BED="final_run_files/bed/expanded_regions_${ENSEMBL_ID}.bed"

# Output files
OUTPUT_VCF="final_run_files/${ENSEMBL_ID}_${BP_DISTANCE}.vcf.bgz"

# Use bcftools to filter VCF by the expanded BED region
chr_vcf=$(bcftools view -H $INPUT_VCF | head -n1 | cut -f1)
chr_bed=$(head -n1 ${EXPANDED_BED} | cut -f1)

# If the VCF chromosome column does not match the format in the interval bed, change it to match.
if [[ "$chr_vcf" != "$chr_bed" ]]; then
    awk -v chr="$chr_vcf" 'BEGIN{OFS="\t"}{$1=chr; print}' "$EXPANDED_BED" > "${EXPANDED_BED}.tmp"
    mv "${EXPANDED_BED}.tmp" "$EXPANDED_BED"
fi

# Note that here, we ensure the correct naming when constructing the gene-specific vcf.
bcftools view --threads "$THREADS" -R "$EXPANDED_BED" "$INPUT_VCF" |
  bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --rename-chrs data/chr_map.tsv |
  bcftools view -Oz -o "${OUTPUT_VCF}"

bcftools index --csi -f "${OUTPUT_VCF}"

echo "Created indexed bgzipped VCF files for the gene ${ENSEMBL_ID}"
