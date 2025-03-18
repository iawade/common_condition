#!/bin/bash

# This will maybe become a .smk

# think about structure
# error handling and loop?
# If there's 70 maybe a snakemake pipeline could work
# Can then collate outputs at the end
# I can't tell if I'm overcomplicating or not
# may want one big conda env yaml or one docker env

# where am I putting all the options? Config yaml?

# 1. Generate gene bed
# Can fail if BioMart server down; Todo error handle

# Run BioMart query to get gene coordinates (can fail if BioMart server down)
## prod 1 line bed - chrom, start, stop
## name of output will be <STABLE_GENE_ID>.bed e.g. ENSG00000186575.bed

python biomart_start_end_query.py --ensembl_id "$ENSEMBL_ID" 

# 2. Identify common variants within specified distance
bash id_variants_for_common_variant_conditioning.sh \
    INPUT_VCF="$1" \
    ENSEMBL_ID="$2" \
    BP_DISTANCE="$3" \
    MAF_COMMON="$4" \
    THREADS


# 3. Create filtered inputs for saige step 2
## can only run one gene at a time as SNPs to condition on change w/ each gene 
grep $ENSEMBL_ID $GROUPFILE > ${ENSEMBL_ID}_saige_group_${trait}.txt
## ok a script would be better as I can add more error handling and keep things neat
## double check that the variants to cindition on don't need to be i nthe group file?

# filter group file

# 4. Run step 2
bash saige_step2_conditioning_check.sh \
    inputs

## Done or echo or whatever

