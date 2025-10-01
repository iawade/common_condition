#!/bin/bash

# Input parameters
INPUT_PLINK="$1" # QC'd, needs GT's; single ancestry-group
ENSEMBL_ID="$2" # Prevents any potential issues with gene symbols
BP_DISTANCE="$3" # TODO Error handling / kb vs just the number
THREADS="$4"
SPARSEGRMID="$5"
OUT_FOLDER="$6"

EXPANDED_BED="${OUT_FOLDER}/bed/expanded_regions_${ENSEMBL_ID}.bed"

# Output files
OUTPUT_PLINK="${OUT_FOLDER}/${ENSEMBL_ID}_${BP_DISTANCE}"

plink2 --bfile ${INPUT_PLINK} \
      --extract bed0 ${EXPANDED_BED} \
      --keep ${SPARSEGRMID} \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 10000 \
      --make-bed \
      --out ${OUTPUT_PLINK} || true

if [ ! -e "${OUTPUT_PLINK}.bim" ]; then
  echo "Edge case - plink file does not exist, following restriction"
  touch ${OUTPUT_PLINK}.bim
  touch ${OUTPUT_PLINK}.bed
  touch ${OUTPUT_PLINK}.fam
else
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
      if ($1 == "23") {
        $1 = "chrX"
      }
      if ($1 == "24") {
        $1 = "chrY"
      }
      print
    }
  }' OFS='\t' ${OUTPUT_PLINK}.bim > ${TMPFILE} && mv ${TMPFILE} ${OUTPUT_PLINK}.bim
fi

echo "Created plink fileset (.bim/.bed/.fam) for the gene ${ENSEMBL_ID}"
