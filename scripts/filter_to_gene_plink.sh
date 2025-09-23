#!/bin/bash

# Input parameters
INPUT_PLINK="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
THREADS="$4"
SPARSEGRMID="$5"

EXPANDED_BED="final_run_files/bed/expanded_regions_${ENSEMBL_ID}.bed"

# First, check to see if this is a superset of the collection of samples
# used to fit the model
n_sparse=$(wc -l < "$SPARSEGRMID")
n_fam=$(wc -l < "${INPUT_PLINK}.fam")

# Output files
OUTPUT_PLINK="final_run_files/${ENSEMBL_ID}_${BP_DISTANCE}"

# Compare and raise error
if (( n_fam < n_sparse )); then
    echo "Error: .fam file ($n_fam samples) has fewer entries than sparse GRM ID file ($n_sparse samples)." >&2
    exit 1
else
    chr_bed=$(head -n1 ${EXPANDED_BED} | cut -f1)
    chr_plink=$(head -n1 ${INPUT_PLINK}.bim | cut -f1)

    # If the VCF chromosome column does not match the format in the bed interval, change the bed to match.
    if [[ "$chr_plink" != "$chr_bed" ]]; then
        awk -v chr="$chr_plink" 'BEGIN{OFS="\t"}{$1=chr; print}' "$EXPANDED_BED" > "${EXPANDED_BED}.tmp"
        mv "${EXPANDED_BED}.tmp" "$EXPANDED_BED"
    fi

    plink2 --bfile ${INPUT_PLINK} \
          --extract bed0 ${EXPANDED_BED} \
          --keep ${SPARSEGRMID} \
          --make-bed \
          --out ${OUTPUT_PLINK}.tmp
    plink2 --bfile ${OUTPUT_PLINK}.tmp \
          --set-all-var-ids @:#:\$r:\$a \
          --make-bed \
          --out ${OUTPUT_PLINK}
    rm ${OUTPUT_PLINK}.tmp.*
fi

TMPFILE=$(mktemp)
awk '{
  if ($1 ~ /^chr/) {
    # already has "chr"
    if ($1 == "chr23") {
      $1 = "chrX"
    }
    if ($1 == "chr24") {
      $1 = "chrY"
    }
    print  
  } else {
    $1 = "chr" $1
    if ($1 == "chr23") {
      $1 = "chrX"
    }
    if ($1 == "chr24") {
      $1 = "chrY"
    }
    print
  }
}' OFS='\t' ${OUTPUT_PLINK}.bim > ${TMPFILE} && mv ${TMPFILE} ${OUTPUT_PLINK}.bim

echo "Created plink fileset (.bim/.bed/.fam) for the gene ${ENSEMBL_ID}"
