#!/bin/bash

# Clear or create the output file first
> all_var_ids.txt

# Loop over all .vcf.gz files in the current directory
for i in *.vcf.gz; do
    [ -e "$i" ] || continue  # Skip loop if no .vcf.gz files exist
    bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$i" >> all_var_ids.txt
done
