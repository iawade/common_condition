#!/bin/bash

# # Activate conda environment
# conda activate biallelic_effects
# source /home/jupyter/anaconda3/etc/profile.d/conda.sh

# TO DO - we're going to make some kind of wrapper that tests a given gene on a given trait (ideally looping over all) 
# I think there's about 70?

## vscode - how to view MD preview as write

## python script for regions
## pull out common variants to condition on script
## wrapper script that performs spa test
## check everything is ready script??
## script where you put in gene, pairs and vcf and get loopin

## check bucket,. etc - use the JSON - also female/male only runs

## coding regions only - prev script

## Need some inputs to test script with 

# needs to be able to work w/ one chromosome or all chrom file
## needs to be able to say - the gene isn't in this file - keep going or move on

## chr1 vs 1 functionality
# check condition format chrom:pos:ref:atl , chr:pos_ref/alt??

## MVP
# think about minimal inputs for step2 

# probably should be testing in UKB...

# filoter group file on ENSG - means spa test will run quick asf

# Input and output variables
VCF="${1}" 
OUT="${2}"
CHR="$3"
MIN_MAC="$4"
MODELFILE="$5"
VARIANCERATIO="$6"
SPARSEGRM="$7"
SPARSEGRMID="${SPARSEGRM}.sampleIDs.txt"
SUBSAMPLES="$8"

step2_SPAtests.R \
        --vcfFile=${VCF} \
        --vcfFileIndex="${VCF}.csi"\
        --vcfField="DS" \
        --chrom="$CHR" \
        --minMAF=0 \
        --minMAC=${MIN_MAC} \
        --GMMATmodelFile=${MODELFILE} \
        --varianceRatioFile=${VARIANCERATIO} \
        --sparseGRMFile=${SPARSEGRM} \
        --sparseGRMSampleIDFile=${SPARSEGRMID} \
        --SampleFile=${SUBSAMPLES} \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --SPAcutoff=0.5 \
        --pCutoffforFirth=0.10 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${OUT} \
        --groupFile=$GROUPFILE \
        --annotation_in_groupTest=$ANNOTATIONS \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --condition="$CONDITION" \
        --maxMissing=1 \
        --impute_method="mean"