# Instructions for confirming BRaVa rare-variant gene-trait associations by conditioning on nearby common variants 

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
      >  - **The pipeline expects files to be named based on the phenotype/trait IDs** from the BRaVa pilot analysis "nominate phenotypes" Google sheet. These IDs are provided in the `phenotype_ID` fields in the JSON. That string just needs to be anywhere in the filenames. For example:

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
      > group_file: "test_files/group_file.txt"
      > ```
> [!NOTE]
>  - Annotations for potential common variants that will be pulled out in this workflow and used for conditioning **are not** required to be in the group file

5. **QCed VCF with Genotypes**
   - VCFs with all samples and sites for analysis included, including individual genotypes
   - Please remember to use a quality-controlled VCF that **still includes common variants**
   - **VCFs must have a `.csi` index**, located in the same place within the file structure as each VCF
   - **Please make sure values/names in VCF `CHROM` column match those in the bed file** (the protein coding intervals BED file, not a plink bed!): i.e. "chr1" instead of "1", "CM000663.2", "NC_000001.11" or something even more esoteric
> [!IMPORTANT]
>  - The pipeline expects a **separate VCF for each chromosome**.

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
##### **plink Conversion**
   - As the SAIGE step2-conditional analysis is on such a small region, the benefits of using plink would be relatively marginal. If variant files are currently in these file formats please convert them first and then use the output VCFs in this workflow
   - Conversion (with bgzip compression and creation of .csi index) can be achieved with `plink1/2` and `bcftools`/`tabix`:
      > #### plink v1.x (.bed/.bim/.fam):
      > ```sh
      > plink --bfile input_prefix --recode vcf --out output_prefix
      > bgzip output_prefix.vcf
      > tabix --csi -p vcf output_prefix.vcf.gz
   
      > #### plink v2 (.pgen/.pvar/.psam):
      > ```sh
      > plink2 --pfile input_prefix --export vcf bgz --out output_prefix
      > bgzip output_prefix.vcf
      > tabix --csi -p vcf output_prefix.vcf.gz

> [!IMPORTANT]
>  - Remember, ensure that the CHROM column in the VCF matches the protein coding BED file format (e.g. "chr1", not "1"). Similarly, ensure the group file also matches ("chr1" not "1") for variants in the group file.

6. **Sparse GRM**
   - From previous analyses there should be a sparse genetic relationship matrix (GRM) generated by SAIGE step0 consisting of the individuals included in the analysis
      > #### `config.yaml`
      > ```yaml
      > sparse_matrix: "test_files/test.sparseGRM.mtx"
      > ```
> [!NOTE]
>  - SPARSE GRM sample IDs are also required for the pipeline. These should be placed in the same location as the sparse GRM file, and must have the following naming convention:
> `{SPARSEGRM_file}.sampleIDs.txt`
> This is the same naming as exported by SAIGE step0, so shouldn't be an issue! 

7. **Other Config Parameters**
   - "gene_trait_pairs_to_test" relates to the csv of BRaVa gene-trait pairs to be tested/confirmed with these conditional analyses. The actual csv is not provided in this repo. It can be downloaded from `gs://brava-meta-pilot-analysis/`. **If you're unable to access as an analyst in BRaVa, please email bravaconsortium@gmail.com.
   - "maf" outlines the minor allele frequencies defining "common" variants to be tested in this conditional analysis
   - "distance" relates to how many bases up and downstream of the gene of interest's start/stop coordinates to identify common variants for conditioning
   - The values stored in these entries **should not be changed** in the config. The config file you use should have:
      > #### `config.yaml` 
      > ```yaml
      > distance: 500000
      > maf:
      >   - 0.001
      >   - 0.0001
      > ```

### Required Software
- **SAIGE**
- **bedtools**
- **bcftools**

#### **Software Handling:**
   - Software can either be independently installed or can all be used via the conda environment provided in this repo 
   
#### Conda Environment Setup
To create the conda environment from the `brava_hits_common_condition_check_conda_env.yaml` file provided, follow these steps:

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

##### **Install bcftools separately
   - bcftools is lightweight and can be installed very easily. See https://www.htslib.org/download/ for instructions.

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

## Notes
- This pipeline utilises `mktemp`. By default this writes files in the `/tmp/` directory. If `/tmp/` is not available, for whatever reason, it may be worth double-checking that `mktemp`'s workarounds (for example writing to `/var/tmp/` ) are working before running the workflow. 
