#!/bin/bash

# Input parameters
INPUT_VCF="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
MAF_COMMON="$4" 
THREADS="$5"

# Output files
VARIANTS_LIST="${ENSEMBL_ID}_${BP_DISTANCE}_${MAF_COMMON}_list.txt"
VARIANTS_COMMA="${ENSEMBL_ID}_${BP_DISTANCE}_${MAF_COMMON}_string.txt"

# Expand the BED regions for query and filter to coding regions
EXPANDED_BED="expanded_regions_${ENSEMBL_ID}.bed"
awk -v BP_DISTANCE="$BP_DISTANCE" 'BEGIN {OFS="\t"} {
    start = $2 - BP_DISTANCE;
    end = $3 + BP_DISTANCE;

    if (start < 0) start = 0;
    
    # Add "chr" prefix if chromosome is numeric
    if ($1 ~ /^[0-9]+$/) $1 = "chr" $1;

    print $1, start, $2;  # Upstream region
    print $1, $3, end;    # Downstream region
}' "$ENSEMBL_ID.bed" | bedtools intersect -a stdin -b protein_coding_regions_hg38_no_padding_no_UTR_v47.bed > "$EXPANDED_BED"

# Use bcftools to filter VCF by the expanded BED regions and MAF threshold
# Using && which is the same as max(MAC > 40, MAF > $MAF_COMMON) ; unless I'm losing the plot
bcftools view --threads "$THREADS" -R "$EXPANDED_BED" "$INPUT_VCF" |
  bcftools filter --threads "$THREADS" -i 'MAC > 40 && MAF > $MAF_COMMON' |
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' >> "$VARIANTS_LIST"

# Use a temporary file for sorted output to avoid in-place modification
SORTED_VARIANTS_LIST=$(mktemp)
sort -V -u "$VARIANTS_LIST" > "$SORTED_VARIANTS_LIST"

# Replace the original list with the sorted list
mv "$SORTED_VARIANTS_LIST" "$VARIANTS_LIST"

# Output for debugging and sanity checking
echo "Variant list written to $VARIANTS_LIST"

# Convert the variants list into a comma-separated format for SAIGE conditioning flag
TEMP_COMMA_FILE=$(mktemp)
awk 'BEGIN {ORS=","} {print $1,$2,$3,$4}' "$VARIANTS_LIST" > "$TEMP_COMMA_FILE"
sed 's/,$/\n/' "$TEMP_COMMA_FILE" > "$VARIANTS_COMMA"

# Cleanup temporary files
rm "$EXPANDED_BED" "$TEMP_COMMA_FILE"

echo "String for --condition flag for SAIGE step 2 written to $VARIANTS_COMMA"
