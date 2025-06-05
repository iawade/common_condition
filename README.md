# Instructions for Confirming BRaVa Gene-Trait "hits" by Conditioning on Nearby Common Variants 
[Please read until the end - some important considerations/potential pitfalls]

## Prerequisites

### Required Files and Config (config.yaml) Set up
1. **JSON File**
   - Provided in this repo (and therefore has no option/entry/key in the config file)
   - It looks like:
      > #### `pilot_phenotypes.json`
      > ```json
      > {
      > "phenotype": "Age-related macular degeneration",
      > "sex_specific_run": "",
      > "phenotype_ID": "AMD",
      > "definition": "ICD_Phecode",
      > "trait_type": "binary",
      > "invnormalise": "FALSE",
      > "tol": 0.02
      > },

2. **Protein-Coding Regions BED File**
   - Any can be used (however please use **Gencode** version **v39** for consistency with the main pilot analysis) or generate one with `bedtools`:
      
      > #### `command-line/sh/zsh`
      > ```sh
      > wget -O - "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz" |\
      > gunzip -c | grep 'transcript_type "protein_coding"' |\
      > awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |\
      > sort -T . -V -t $'\t' -k1,1 -k2,2n | bedtools merge > protein_coding_regions_hg38_no_padding_no_UTR_v39.bed
      > ```

   - Which should look like:
      > #### `protein_coding_regions_hg38_no_padding_no_UTR_v39.bed`
      > ```
      > chr1    65418   65433
      > chr1    65519   65573
      > chr1    69036   71585
      > chr1    450739  451678
      > chr1    685715  686654
      > ```
   - Note: if you don't have `bedtools` already installed it can be installed as part of the conda environment provided in this repo (see below)

3. **Model and Variance Ratio Files**
   - GMMAT model files (are given `.rda` extension by default in SAIGE step1) from prior BRaVa pilot analysis
   - Variance ratio files (given `varianceRatio.txt` extension in SAIGE step1)
   - The config file (config.yaml) has keys/entries for a list of model files and variance ratio files to be used, locations relative to where the pipeline is being run from - change input name as needed
      > #### `config.yaml`
      > ```yaml
      > list_of_model_files: "model_file_list.txt"
      > variance_ratio_file: "variance_file_list.txt"
      > ```
   - **The pipeline expects files to be named based on the phenotype/trait ID's** from the BRaVa pilot analysis "nominate phenotypes" Google sheet. These ID's are provided in the `phenotype_ID` fields in the JSON. That string just needs to be anywhere in the filenames; for example:
      > #### `variance_file_list.txt`
      > ```txt
      > test_files/TG.varianceRatio.txt
      > test_files/Urolith.varianceRatio.txt
      > test_files/but_they_can_be_AMD_called_anything.extension
      > ```

      > #### `model_file_list.txt`
      > ```txt
      > test_files/TG.rda  
      > test_files/Urolith.rda
      > ```
       
4. **QC’d Group File (see warning below)**
   - Contains all annotated variants used for the pilot analysis
   - The pipeline currently only supports one group file - if previously split by chromosome please concatenate together
   - Two lines per gene, one with a space-delineated list of variants, followed by the corresponding annotations:
      > #### `group file (name and extension are unimportant)`
      > ```
      > ENSG00000172967 var chr22:16783569:T:A chr22:16783577:A:G chr22:16783578:G:A chr22:16783579:C:T 
      > ENSG00000172967 anno non_coding synonymous pLoF damaging_missense
      > ```
   
      > #### `config.yaml`
      > ```yaml
      > group_file: "test_files/group_file.txt"
      > ```
   - Annotations for potential common variants that will be pulled out in this workflow and used for conditioning **are not** required to be in the group file


5. **QC’d VCF with Genotypes**
   - VCFs with all samples and sites for analysis included, including individual genotypes
   - Please remember to use a quality-controlled VCF that **still includes common variants**
   - **VCFs must have a `.csi` index**, located in the same place within the file structure as the/each VCF
   - **Please make sure values/names in VCF `CHROM` column match those in the bed file**: i.e. "chr1" instead of "1", "CM000663.2", "NC_000001.11" or something even more esoteric
   - Multiple (although untested as to whether it will work if VCFs are chunked into multiple VCFs per chromosome) or a single VCF are supported. Either option requires a file with the location of VCF[s]
      > #### `config.yaml`
      > ```yaml
      > list_of_vcf_files: "vcf_list.txt"
      > ```
      
      > #### `vcf_list.txt`
      > ```
      > test_files/QC_applied_common_variants_present_chr21.vcf.bgz
      > or
      > ...
      > test_files/QC_applied_common_variants_present_chr1.vcf.bgz
      > test_files/QC_applied_common_variants_present_chr2.vcf.bgz
      > ```
---
##### **Converting PLINK or BGEN to VCF**

This workflow requires input variant files in VCF format. If your data is currently in PLINK (`.bed/.bim/.fam` or `.pgen/.pvar/.psam`) or BGEN format, **you must convert it to VCF first**.

   ##### Important Notes on PLINK 1.x

   PLINK v1.9 can alter data in ways that cause serious issues:

   - Chromosome names may be output as numeric (e.g. `1`) instead of `chr1`
   - Sample IDs in the VCF may include family IDs (e.g. `FID_IID`)
   - Reference/alternate allele order may be flipped unless `--keep-allele-order` is used

   An example script is provided to safely convert from PLINK v1 format:
      scripts/UTILITY_export_vcf_from_plink.sh

      This script:
      - Enforces `chr1`, `chr2`, ... chromosome naming
      - Strips family IDs from sample names
      - Uses `--keep-allele-order` to preserve allele direction
      - Outputs a bgzipped VCF with `.csi` index
---


6. **Sparse GRM**
   - From previous analyses there should be a sparse genetic relationship matrix (GRM) generated by SAIGE step0 consisting of the individuals included in the analysis
      > #### `config.yaml`
      > ```yaml
      > sparse_matrix: "test_files/test.sparseGRM.mtx"
      > ```

7. **Other Config Parameters**
   - "gene_trait_pairs_to_test" relates to the csv of BRaVa gene-trait pairs to be tested/confirmed with these conditional analyses. The actual csv is not provided in this repo. It can be downloaded from `gs://brava-meta-pilot-analysis/`
   - "maf" outlines the minor allele frequencies defining "common" variants to be tested in this conditional analysis
   - "distance" relates to how many bases up and downstream of the gene of interest's start/stop coordinates to identify common variants for conditioning
   - The values stored in these entries do not need to be changed in the config 

### Required Software
- **SAIGE**
- **bedtools**
- **BioMart for Python**
- **bcftools**

#### **Software Handling:**
   - Software can either be independently installed or can all be used via the conda environment provided in this repo 
   
#### Conda Environment Setup
To create the Conda environment from the `brava_hits_common_condition_check_conda_env.yaml` file provided, follow these steps:

##### **Activate Conda:**
   - Hopefully conda is already installed but, if necessary, it can be installed with instructions from https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
   - Activate Conda (double-check with system admin if unsure how)
   - For example it can activated with:
      > ```bash
      > eval "$(conda shell.bash hook)"
      > ```
   - Or on some systems:
      > ```bash
      > "/opt/software/applications/anaconda/3/etc/profile.d/conda.sh"
      > ```

##### **Create the Conda Environment:**
   - Run the following command to install the environment in the desired location:
     > ```bash
     > conda env create --prefix /path/to/envs/brava_hits_common_condition_check -f brava_hits_common_condition_check_conda_env.yaml
     > ```
   - Note: This is a relatively large environment and needs and a machine with at least 4-8Gb of memory to 

##### **Activate the Environment:**
   - Once the environment is created, activate it using the following command:
     > ```bash
     > conda activate /path/to/envs/brava_hits_common_condition_check
     > ```

##### **Note on bcftools and libcrypto (especially if using provided Conda environment):**
  
- **libcrypto Symlink for Compatibility**: If there are issues with `libcrypto` and `bcftools` compatibility, a symlink can be created to bypass package conflicts when using bcftools with conda. Here's how the symlink can be created:

  1. Locate the `libcrypto` library in the Conda environment:
     > ```bash
     > find /path/to/envs/brava_hits_common_condition_check/ -name "*libcrypto*"
     > ```
  2. Which will output something like:
     > ```
     > /path/to/envs/brava_hits_common_condition_check/lib/libcrypto.so
     > /path/to/envs/brava_hits_common_condition_check/lib/libcrypto.so.3
     > ```
  3. Create a symlink from `libcrypto.so.1.0.0` to the version available in the environment:
     > ```bash
     > ln -s /path/to/envs/brava_hits_common_condition_check/lib/libcrypto.so.3 /path/to/envs/brava_hits_common_condition_check/lib/libcrypto.so.1.0.0
     > ```
  4. Verify the symlink:
     > ```bash
     > ls -l /path/to/envs/brava_hits_common_condition_check/lib/libcrypto.so.1.0.0
     > ```
   Note:
   - TODO 
   - if still not working need a different symlink - rm prev link
    -- ldconfig -p | grep libcrypto
        libcrypto.so.1.1 (libc6,x86-64) => /lib/x86_64-linux-gnu/libcrypto.so.1.1
        libcrypto.so (libc6,x86-64) => /lib/x86_64-linux-gnu/libcrypto.so
   -- do rm envs/brava_hits_common_condition_check/lib/libcrypto.so.1.0.0 
 ln -s /lib/x86_64-linux-gnu/libcrypto.so.1.1 envs/brava_hits_common_condition_check/lib/libcrypto.so.1.0.0

  This step ensures that the required version of `libcrypto` is linked correctly and avoids issues with package compatibility.

- **Alternative**: If `bcftools` is already installed , it can simply be loaded (for example with `module load BCFtools`) to run the pipeline, bypassing the need to create the symlink.

## Running the SnakeMake pipeline
- Once everything is set up, running should be relatively straightforward
- A submission script `snakemake_iterative_call.sh` is provided, simply run the script locally or submit to a cluster with relevant options, for example:
   > ```bash
   > bash snakemake_iterative_call.sh
   > ```
   or
   > ```bash
   > sbatch --job-name=common_conditioning_pipeline --mem-per-cpu=8000  --ntasks=1 --cpus-per-task=8 --partition=short --output=slurm-%x-%A_%a.out snakemake_iterative_call.sh
   > ```
   etc

## Final Outputs
- The pipeline ultimately produces two files:
   - `brava_stepwise_conditional_analysis_results.txt`
   - `brava_stepwise_conditional_analysis_results.txt.singleAssoc.txt`
- These contain all the conditional SAIGE step2 outputs concatenated together
- Please rename to something appropriate and send to Duncan Palmer!

---

## Group File Reordering — CRITICAL WARNING

   WARNING: Variant order must exactly match between your VCF and the group file.  
   If not, SAIGE will either silently skip variants or misalign them, leading to incorrect or silently corrupted results.

   SAIGE reads variants from the group file in the order they appear, assuming the VCF files contain variants in the exact same order. If the order does not match:

   - Variants may be skipped or misread
   - SAIGE may assign the wrong variant annotations
   - Your association results could be incomplete or wrong, with no explicit error message

   This is especially important if:
   - You've subsetted or filtered the VCF
   - You merged or re-sorted group files
   - You're running with partial data from the full variant list

   ### Recommended Fix

   We provide utility scripts to regenerate a reordered and filtered group file to match the actual VCF content.

   Step 1: Extract variant IDs from your working VCF(s)
   Step 1: Extract variant IDs from your working VCF(s)

      bash scripts/UTILITY_extract_chrom_pos_ref_alt_from_vcf.sh

   This creates a file `all_var_ids.txt` with one line per variant in the format:

      CHROM:POS:REF:ALT

   Step 2: Reorder the group file based on those VCF variants

      python scripts/UTILITY_reorder_group_file_based_on_extracted_vars_from_vcf.py

   This script:
   - Filters out variants not present in the VCF
   - Keeps the original order from the VCF
   - Writes a new `filtered_group_file.txt.gz`

   ### Expected Inputs

   - `groups/combined_group_file.txt.gz` — your original group file
   - `*.vcf.gz` — the VCFs you're using for association tests

   If your paths differ, edit the paths in the scripts accordingly.

   - The script assumes variants are listed as `CHROM:POS:REF:ALT` in both the VCF and group file.
   - If this format differs in your data, you'll need to modify the extraction and matching logic accordingly.
   - Keep the `var` and `anno` lines in sync — they must remain paired after filtering.

   ### Final Reminder

   SAIGE will not warn you if variant IDs mismatch or go out of order.  
   You must manually ensure correct ordering, or you risk invalid results.
   
   ---

## Notes
- This pipeline utilises `mktemp`. By default this writes files in the `/tmp/` directory. If `/tmp/` is not available, for whatever reason, it may be worth double-checking that `mktemp`'s workarounds (for example writing to `/var/tmp/` ) are working before running the workflow 