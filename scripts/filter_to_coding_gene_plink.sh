#!/bin/bash

# Input parameters
INPUT_PLINK="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
MAF_COMMON="$4" 
THREADS="$5"
SPARSEGRMID="$6"

EXPANDED_BED="run_files/bed/expanded_regions_${ENSEMBL_ID}.bed"

# First, check to see if this is a superset of the collection of samples
# used to fit the model
n_sparse=$(wc -l < "$SPARSEGRMID")
n_fam=$(wc -l < "${INPUT_PLINK}.fam")

# Output files
OUTPUT_PLINK="run_files/${ENSEMBL_ID}_${BP_DISTANCE}_${MAF_COMMON}"

# Compare and raise error
if (( n_fam < n_sparse )); then
    echo "Error: .fam file ($n_fam samples) has fewer entries than sparse GRM ID file ($n_sparse samples)." >&2
    exit 1
else
    # Use plink to filter by the expanded BED regions and MAF threshold
    # Using --maf $MAF_COMMON and --mac 41, which is the same as max(MAC > 40, MAF > $MAF_COMMON)
    echo "Filtering PLINK files (${INPUT_PLINK}.*) to match sparse GRM IDs ($n_sparse samples)..."
    plink2 --bfile ${INPUT_PLINK} \
      --extract bed0 ${EXPANDED_BED} \
      --keep ${SPARSEGRMID} \
      --maf ${MAF_COMMON} \
      --mac 41 \
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
    print  
  } else {
    $1 = "chr" $1
    if ($1 == "chr23") {
      $1 = "chrX"
    }
    print
  }
}' OFS='\t' ${OUTPUT_PLINK}.bim > ${TMPFILE} && mv ${TMPFILE} ${OUTPUT_PLINK}.bim

echo "Created plink fileset (.bim/.bed/.fam) for the gene ${ENSEMBL_ID}"
