# Instructions for confirming BRaVa rare-variant gene-trait associations by conditioning on nearby common variants 
> [!IMPORTANT]
> Please read until the end - some important considerations/potential pitfalls

## Motivation and outline of the procedure

This pipeline is designed to determine common variants which have the potential to be driving rare variant associations uncovered in our meta-analysis. We will then condition on these variants in a final collection of gene-based association tests at candidate (gene, phenotype) pairs.

Briefly, following gene-based analysis and subsequent meta-analysis, we will carry out common variant conditioning for all gene-phenotype pairs with significant meta-analysis association _P_-values. For each (biobank, ancestry, phenotype) tuple, we will then determine a collection of variants to condition on to account for signal driven by common variation nearby. To do this, we have created the iterative conditioning pipeline in this repository (`snakemake_iterative_call.sh`).

For each MAF mask, we carry out association analysis of all variants with MAF greater than the MAF of the mask within 500kb of the gene. If any variant association has an association _P_-value < 1 × 10<sup>-5</sup>, we add it to a set of conditioning variants and condition on it `--condition` flag within SAIGE, and iteratively rerun until no variant in the region is associated (_P_-value < 1 × 10<sup>-5</sup>) with the trait. This procedure is carried out for all (ancestry, biobank) pairs.

Each biobank (you!) then provides conditioning variant lists, as well as gene and variant based association test results for subsequent testing of agreement with earlier summary statistic results.

> [!WARNING]  
> In the config, make sure that for `annotations_to_include` that the damaging missense and other missense naming is as in your annotation group file. This is likely one of `damaging_missense_or_protein_altering` or `damaging_missense`. Similarly for other missense: likely one of `other_missense_or_protein_altering` or `other_missense`. Multiple naming conventions were used in the initial return of sumstats! The current default in the example config in this directory is
> ```
> annotations_to_include: "pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous"
> ```
> Please check!

In a final step, we then determine the union of these lists centrally for each genetic ancestry. The resultant variant lists are then shared back with the constituent biobanks. Biobanks then perform final gene-based association analysis conditioning on these variants. To guard against collinearity in the variants used for conditioning, we first perform linkage disequilibrium pruning, ensuring that no-pair of variants in the set have _r_<sup>2</sup> > 0.9. (WORK IN PROGRESS).

## Prerequisites

### Required Files and Config (config.yaml - see configs folder for examples) set up
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
   - Provided within this repo in the data folder. Generated using **Gencode** version **v39** for consistency with the main pilot analysis. Code used to generate the file is:
      
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
   - The config files (.yaml files within the configs folder) have example keys/entries for a list of model files and variance ratio files to be used, locations are relative to where the pipeline is being run from - change input name as needed
      > #### `config.yaml`
      > ```yaml
      > list_of_model_files: "model_file_list.txt"
      > variance_ratio_file: "variance_file_list.txt"
      > ```
> [!IMPORTANT]
> **The pipeline expects files to be named based on the phenotype/trait IDs** from the BRaVa pilot analysis "nominate phenotypes" Google sheet. These IDs are provided in the `phenotype_ID` fields in the JSON. That string just needs to be anywhere in the filenames.
   - For example:
      > #### `variance_file_list.txt`
      > ```txt
      > test_files/TG.varianceRatio.txt
      > test_files/Urolith.varianceRatio.txt
      > test_files/but_they_can_be_AMD_called_anything_but_must_end_with.varianceRatio.txt
      > ```
      > #### `model_file_list.txt`
      > ```txt
      > test_files/TG.rda  
      > test_files/Urolith.rda
      > ```
       
4. **QCed Group File**
   - Contains all annotated variants used for the pilot analysis (this should be the group file used in the pilot analysis)
   - The pipeline currently only supports one group file - if previously split by chromosome please concatenate together
   - As a reminder, there are two lines per gene, one with a space-delineated list of variants, followed by the corresponding annotations:
      > #### `group file (name and extension are unimportant)`
      > ```
      > ENSG00000172967 var chr22:16783569:T:A chr22:16783577:A:G chr22:16783578:G:A chr22:16783579:C:T 
      > ENSG00000172967 anno non_coding synonymous pLoF damaging_missense
      > ```
   
      > #### `config.yaml`
      > ```yaml
      > list_of_group_files: "test_files/group_file_list.txt"
      > ```
> [!NOTE]
>  - Annotations for potential common variants that will be pulled out in this workflow and used for conditioning **are not** required to be in the group file

5. **QCed VCF or plink file with Genotypes**
   - plink (`.bim`/`.bed`/`.fam`) files or VCF files with all samples and sites for analysis included, including individual genotypes
   - Please remember to use a quality-controlled VCF/plink files that **still includes common variants**
   - **VCFs must have a `.csi` index**, located in the same place within the file structure as each VCF
   - **Please make sure values/names in VCF `CHROM` column match those in the bed file** (the protein coding intervals BED file, not a plink bed!): i.e. "chr1" instead of "1", "CM000663.2", "NC_000001.11" or something even more esoteric.
   - Note that we actually try to deal with the above, but things will work faster if they're in the same format
> [!IMPORTANT]
> The pipeline expects a **separate genetic data file (plink or vcf) for each chromosome**.
> Additionally, these **must be named with a "." either side of the chrom name** 
   - For example, in the case of vcf files being passed:
     #### `config.yaml`
     ```yaml
     list_of_vcf_files: "vcf_list.txt"
     ```
     #### `vcf_list.txt`
     ```
     test_files/QC_applied_common_variants_present.chr1.vcf.bgz
     test_files/QC_applied_common_variants_present.chr2.vcf.bgz
     ```
   - Similarly, if plink files are passed:
     #### `config.yaml`
     ```yaml
     list_of_plink_files: "plink_list.txt"
     ```
     #### `plink_list.txt`
     ```
     test_files/QC_applied_common_variants_present.chr1
     test_files/QC_applied_common_variants_present.chr2
     ```
Note that in the case if plink files, the extension is excluded, but it is expected that the `.bim`, `.bed`, and `'.fam` are in the same location, with the same naming before the extension.

---
#### **Converting to VCF or plink (`.bim`/`.bed`/`.fam`)**

This workflow requires input variant files in VCF format or plink (`.bed/.bim/.fam`) format. If your data is currently in plink2 (`.pgen/.pvar/.psam`) or BGEN format, **you will need to convert it to plink or VCF first** - we recommend plink, given the choice, as the SAIGE codebase is more stable for plink files.

##### Important Notes on PLINK 1.x
> [!WARNING]  
> PLINK v1.9 can alter data in ways that cause issues:

   - Chromosome names may be output as numeric (e.g. `1`) instead of `chr1`
   - Sample IDs in the VCF may include family IDs (e.g. `FID_IID`)
   - Reference/alternate allele order may be flipped unless `--keep-allele-order` is used
---

> [!IMPORTANT]
> Remember, ensure that the CHROM column in the VCF or plink files match the protein coding BED file format (e.g. "chr1", not "1"). Similarly, ensure the group file also matches ("chr1" not "1") for variants in the group file. Note that we do our best to check and fix this along the way, but despite my best efforts to correct all edge cases and munge variant IDs to the expected formats, errors could slip through. If they do, it will likely be at the final group-based testing phase, and you'll see the error manifest as a difference in p-value for the unconditioned gene-based test compared to your first run through with universal-saige. I will perform this sanity check, and we can work through any issues if/when they arise.

6. **Sparse GRM**
   - From previous analyses there should be a sparse genetic relationship matrix (GRM) generated by SAIGE step0 consisting of the individuals included in the analysis
      > #### `config.yaml`
      > ```yaml
      > sparse_matrix: "test_files/test.sparseGRM.mtx"
      > ```
> [!NOTE]
> SPARSE GRM sample IDs are also required for the pipeline. These should be placed in the same location as the sparse GRM file, and must have the following naming convention:
> `{SPARSEGRM_file}.sampleIDs.txt`
> This is the same naming as exported by SAIGE step0, so shouldn't be an issue! 

7. **Other Config Parameters**
   - "input_format" should be either `plink` or `vcf`.
   - "gene_trait_pairs_to_test" relates to the csv of BRaVa gene-trait pairs to be tested/confirmed with these conditional analyses. **The actual csv is not provided in this repo**. It can be downloaded from `gs://brava-meta-pilot-analysis/gene_phenotype_pairs_101025.csv.gz`. *If you're unable to access as an analyst in BRaVa, please email bravaconsortium@gmail.com*.
   - "maf" outlines the minor allele frequencies defining "common" variants to be tested in this conditional analysis
   - "distance" relates to how many bases up and downstream of the gene of interest's start/stop coordinates to identify common variants for conditioning
   - The values stored in these entries **should not be changed** in the config. The config file you use should have:
      > #### `config.yaml` 
      > ```yaml
      > conditioning_pvalue: 0.00001
      > distance: 500000
      > maf:
      >   - 0.001
      >   - 0.0001
      > ```
---

## Container Setup (Docker or Singularity/Apptainer)
To create the everything required to run the iterative conditioning pipeline, first determine whether you are able to use docker, or singuarity/apptainer. Once you know which you can use, move to the relevant section below:

### Apptainer (previously known as apptainer)
You can grab the container from dockerhub or google artifact registry. They are identical.
#### From dockerhub
```sh
singularity pull brava-common-check.sif docker://astheeggeggs/brava-common-check:v4
# or
apptainer pull brava-common-check.sif docker://astheeggeggs/brava-common-check:v4
```
#### From google artifact registry
```sh
singularity pull brava-common-check.sif docker://gcr.io/weighty-elf-452116-c7/brava-common-check:v4
# or
apptainer pull brava-common-check.sif docker://gcr.io/weighty-elf-452116-c7/brava-common-check:v4
```
Note that there's a preference for dockerhub (as it's free!), but gcr was more stable in our hands when using dsub on the All of Us platform, so we have provided the two locations.

If you have errors to do with disk quotas - it could well be the your local home directory is very small. You can change the location of the cache to somewhere on your HPC with more available space using e.g.
```sh
mkdir -p .apptainer/cache
export SINGULARITY_CACHEDIR=.apptainer/cache
export APPTAINER_CACHEDIR=.apptainer/cache
```
where you can replace `.apptainer/cache` with where you would like to place it.

### Docker
#### From dockerhub
```sh
docker pull astheeggeggs/brava-common-check:v4
```
#### From google artifact repository
```sh
docker pull gcr.io/weighty-elf-452116-c7/brava-common-check:v4
```

> [!WARNING]
> Depending on your system’s connectivity or network speed, building the Docker/Singularity image can take over an hour.

## Running the pipeline!
Once you have the docker or singularity container pulled, you can then run step 1 of the pipeline. First, you need to clone this repository e.g.:
```
git clone https://github.com/astheeggeggs/common_condition.git
```
You can then run `run_step1_iterative_conditioning.sh` from within the docker that you've pulled. If you're using docker:
```sh
docker run \
      -v $(pwd)/common_condition:/common_condition \
      -w /common_condition \
      gcr.io/weighty-elf-452116-c7/brava-common-check:v4 \
      bash /common_condition/run_step1_iterative_conditioning.sh
```
or if you're using singularity/apptainer:
```sh
apptainer exec \
  -B $(pwd)/common_condition:/common_condition \
  -W /common_condition \
  brava-common-check.sif \
  bash /common_condition/run_step1_iterative_conditioning.sh
```
or
```
singularity exec \
  -B $(pwd)/common_condition:/common_condition \
  -W /common_condition \
  brava-common-check.sif \
  bash /common_condition/run_step1_iterative_conditioning.sh
```

Importantly, all of your data (paths will be in the `config.yml`, `list_of_plink_files`/`list_of_vcf_files`, `list_of_model_files`, `list_of_variance_ratio_files`, `sparse_matrix`, `list_of_group_files`, `gene_trait_pairs_to_test`, `protein_coding_region_bed`, and `phenotype_json`. Each location must be present within the mounted directory which is passed to the container, and filepaths will be relative to the working directory (`/common_condition` in the examples above). I've found that a straightforward approach is to place a folder containing your data this directory and mount it (likely named common_condition), by `-v`/`-w` or `-B`/`-W` as above.

## Final Outputs
- The pipeline ultimately produces three files:
   - `brava_stepwise_conditional_analysis_results.txt`
   - `brava_stepwise_conditional_analysis_results.txt.singleAssoc.txt`
   - `brava_stepwise_conditional_analysis_results.txt.conditioning.variants.txt`
- These contain all the conditional SAIGE step2 outputs and the gene and variant level, as well as the collection of conditioning variants used for each (gene, trait, max_MAF) tuple (where conditioning was applied).
- Please rename to something appropriate and send to Duncan Palmer!

---

### Group File Reordering if using VCF files

> [!WARNING]  
> If you are running the pipeline using VCF files, variant order must exactly match between your VCF and the group file.  
> If not, SAIGE will either silently skip variants or misalign them, leading to incorrect or silently corrupted results.

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

   #### Expected Inputs

   - `groups/combined_group_file.txt.gz` — your original group file
   - `*.vcf.gz` — the VCFs you're using for association tests

   If your paths differ, edit the paths in the scripts accordingly.

   - The script assumes variants are listed as `CHROM:POS:REF:ALT` in both the VCF and group file.
   - If this format differs in your data, you'll need to modify the extraction and matching logic accordingly.
   - Keep the `var` and `anno` lines in sync — they must remain paired after filtering.

   #### Final Reminder

   SAIGE will not warn you if variant IDs mismatch or go out of order.  
   You must manually ensure correct ordering, or you risk invalid results.
   
   ---
