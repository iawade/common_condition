#!/bin/bash

# TODO Some intro text

# just want to check this works - then will make it so input is gene+trait pair?

# Input parameters
INPUT_VCF="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
MAF_COMMON="$4" 

#  TODO make dynamic
# Output files
VARIANTS_LIST="variants_list.txt"
VARIANTS_COMMA="variants_comma.txt"

# TODO this is the only bit that needs an env for people; Can fail if BioMart server down; Todo error han

eval "$(conda shell.bash hook)"
conda activate biomart # conda env for python/biomart script/ handling manual start and end coords?

# Step 1: Run BioMart query to get gene coordinates (can fail if BioMart server down)
python biomart_start_end_query.py --ensembl_id "$ENSEMBL_ID" # prod 1 line bed - chrom, start, stop

# bcftools get variants up and downstream of coords by x bp pre and after
# Then filter MAF 

# Step 2: Expand the BED regions for query
EXPANDED_BED="expanded_regions.bed"
awk -v BP_DISTANCE="$BP_DISTANCE" 'BEGIN {OFS="\t"} {
    start = $2 - BP_DISTANCE;
    end = $3 + BP_DISTANCE;

    if (start < 0) start = 0;
    
    # Add "chr" prefix if chromosome is numeric
    if ($1 ~ /^[0-9]+$/) $1 = "chr" $1;

    print $1, start, $2;  # Upstream region
    print $1, $3, end;    # Downstream region
}' "$ENSEMBL_ID.bed" > "$EXPANDED_BED"

# TODO add multithreadedness
# Step 3: Use bcftools to filter VCF by the expanded BED regions and MAF threshold
bcftools view -R "$EXPANDED_BED" "$INPUT_VCF" |
  bcftools filter -i "MAF>$MAF_COMMON" |
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > "$VARIANTS_LIST"
# Output for debugging and sanity checking

echo "Variant list written to $VARIANTS_LIST"

# Step 4: Convert the variants list into a comma-separated format for SAIGE conditioning flag
awk 'BEGIN {OFS=":"} {print $1,$2,$3,$4}' "$VARIANTS_LIST" |
  tr '\n' ',' | sed 's/,$/\n/' > "$VARIANTS_COMMA"

# Cleanup temporary file
# rm "$EXPANDED_BED" - uncomment again just testing

echo "String for --condition flag for SAIGE step 2 written to $VARIANTS_COMMA"
