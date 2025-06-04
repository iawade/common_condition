#!/bin/bash

# Usage:
#   bash UTILITY_export_vcf_from_plink.sh <PLINK_PREFIX> [THREADS]
#
# Example:
#   bash UTILITY_export_vcf_from_plink.sh UKB.chr3.full_qc.exome.EUR 32
#
# This script converts a PLINK file to a bgzipped, indexed VCF.
# It assumes:
#   - The input PLINK prefix points to .bed/.bim/.fam files.
#   - Chromosomes in the plink files are numeric (e.g., 1..22 or X), without "chr" prefix.
#   - Family IDs in .fam are all zero, resulting in sample IDs like "0_12345" — the "0_" is removed.
#     If FIDs are nonzero, you can change the `gsub(/^0_/, "", $i)` to `gsub(/^[^_]*_/, "", $i)`
#     to remove any FID.

set -euo pipefail

# --- Args
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <PLINK_PREFIX> [THREADS]"
  exit 1
fi

plink_prefix="$1"
threads="${2:-8}"

# --- Derive names
out_prefix="${plink_prefix}.convert"
tmp_vcf="${out_prefix}.tmp.vcf"

# Extract CHR for sed (assumes format like chr3/full_qc.exome.*)
chr=$(basename "$plink_prefix" | sed -E 's/.*chr([0-9XY]+).*/\1/')

echo "[INFO] Running with:"
echo "       PLINK file prefix : $plink_prefix"
echo "       Chromosome        : $chr"
echo "       Threads           : $threads"

# --- Convert to VCF
plink --bfile "${plink_prefix}" \
      --recode vcf \
      --keep-allele-order \
      --out "${tmp_vcf%.vcf}"

# --- Add chr prefix and clean sample IDs
sed "s/^${chr}\b/chr${chr}/" "${tmp_vcf}" > "${tmp_vcf}.chr"

awk 'BEGIN{OFS="\t"} /^#CHROM/{for(i=10;i<=NF;i++){gsub(/^0_/,"",$i)}; print; next} {print}' "${tmp_vcf}.chr" \
  | bgzip -@ "${threads}" -c > "${out_prefix}.vcf.gz"

# --- Index
bcftools index --csi --threads "${threads}" "${out_prefix}.vcf.gz"

# --- Cleanup
rm -f *.log *.nosex "${tmp_vcf}" "${tmp_vcf}.chr"

echo "[DONE] VCF written: ${out_prefix}.vcf.gz"
