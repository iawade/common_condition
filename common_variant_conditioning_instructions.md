# Instructions for Confirming BRaVa Gene-Trait "hits" by Conditioning on Nearby Common Variants 

## Prerequisites

### Required Files
- **JSON File**: Contains model and variance data for each trait. Each trait's files must be named according to the `phenotype_ID` field in the provided JSON.
- **Protein-Coding Regions BED File**: Use your own or generate one with `bedtools`:
  ```sh
  wget -O - "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz" |\
  gunzip -c | grep 'transcript_type "protein_coding"' |\
  awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |\
  sort -T . -V -t $'\t' -k1,1 -k2,2n | bedtools merge > protein_coding_regions_hg38_no_padding_no_UTR_v47.bed
  ```
- **Model and Variance Ratio Files**: Previously generated and named with `phenotype_ID`, stored in the main directory.
- **QC’d Group File**: Contains annotated variants from prior analyses.
- **QC’d VCF with Genotypes**: The same VCF used in prior analysis, with a `.csi` index.

### Required Software
- **SAIGE**
- **bedtools**
- **biomart for Python**

### Conda Environment Setup

To create the Conda environment from the `brava_hits_common_condition_check_conda_env.yaml` file, follow these steps:

### **Create the Conda Environment:**
   - Run the following command to install the environment in your desired location:
     ```bash
     conda env create --prefix /path/to/your/envs/brava_hits_common_condition_check -f brava_hits_common_condition_check_conda_env.yaml
     ```

### **Activate the Environment:**
   - Once the environment is created, activate it using the following command:
     ```bash
     conda activate /path/to/your/envs/brava_hits_common_condition_check
     ```