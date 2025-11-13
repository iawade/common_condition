# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis
# Supports both VCF and PLINK input formats

configfile: "config.yaml"

wildcard_constraints:
    gene="[A-Za-z0-9]+",
    trait="[A-Za-z0-9]+",
    gene_trait="[A-Za-z0-9]+_[A-Za-z0-9]+",  # For concatenated gene_trait wildcards
    distance="[0-9]+",
    maf="[0-9.]+"

# # Testing
# input_format = "plink"
# list_of_input_files = "ukb_paths/eur/eur_plink_file_list.txt"
# list_of_loco_files = "ukb_paths/eur/eur_loco_file_list.txt"
# list_of_regenie_annotation_files = "ukb_paths/eur/eur_annotation_file_list.txt"
# list_of_regenie_setlist_files = "ukb_paths/eur/eur_setlist_file_list.txt"
# keep = "snakemake_eur/keep/ukb22418_b0_v2.autosomes.qced.EUR.id"
# phenotype_file = "snakemake_eur/phenotypes/ukb.standing_height.20250508.tsv"
# covariate_file = "snakemake_eur/covariates/ukb_brava_default_covariates.20250508.tsv"
# gene_trait_pairs_to_test = "gene_phenotype_pairs_051125.csv"
# protein_coding_region_bed = "data/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed"
# phenotype_json = "pilot_phenotypes.json"
# covariate_cols = "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
# categorical_covariate_cols = "sex"
# gene_trait_pairs_to_test = "gene_phenotype_pairs_051125.csv"
# protein_coding_region_bed = "data/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed"
# phenotype_json = "pilot_phenotypes.json"
# distance=500000
# maf=0.001
# min_mac=0.5
# annotations_to_include="pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous"

# # Additional files
# phenotypes="snakemake_eur/phenotypes/ukb.standing_height.20250508.tsv"
# covariates="snakemake_eur/covariates/ukb_brava_default_covariates.20250508.tsv"

# # End testing

# gene_trait_pairs_to_test: "gene_phenotype_pairs_051125.csv"
# protein_coding_region_bed: "data/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed"
# phenotype_json: "pilot_phenotypes.json"

# Read inputs from config file
input_format = config["input_format"] # Has to be plink
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
phenotype_json = config["phenotype_json"]

# Phenotypes should be defined as the collection of available columns in the phenotype file
# intersect that with the json to define the collection of available phenotypes
phenotype_file = config["phenotype_file"]
covariate_file = config["covariate_file"]
keep = config["keep"]

# Format-specific file lists
if input_format == "plink":
    list_of_input_files = config["list_of_plink_files"]
else:
    raise ValueError(f"Invalid input_format: {input_format}. Must be 'plink'")

list_of_loco_files = config["list_of_loco_files"]
list_of_regenie_annotation_files = config["list_of_regenie_annotation_files"]
list_of_regenie_setlist_files = config["list_of_regenie_setlist_files"]

# Load input files
with open(list_of_input_files) as f:
    input_files = [line.strip() for line in f]

# Create format-specific file lists for PLINK
if input_format == "plink":
    plink_bim_files = [f"{p}.bim" for p in input_files]
    plink_bed_files = [f"{p}.bed" for p in input_files]
    plink_fam_files = [f"{p}.fam" for p in input_files]

# Load group files
with open(list_of_regenie_annotation_files) as f:
    regenie_annotation_files = [line.strip() for line in f]

with open(list_of_regenie_setlist_files) as f:
    regenie_setlist_files = [line.strip() for line in f]

print(f"Input format: {input_format}")
print(f"Gene trait pairs to test: {gene_trait_pairs_to_test}")
print(f"Input files Sample: {input_files[:5]}")
print(f"Group files Sample: {regenie_annotation_files[:5]}")
print(f"Setlist files Sample: {regenie_setlist_files[:5]}")

distance = config["distance"]
maf = config["maf"]
min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re
from scripts.extract_chromosome import get_gene_chr
from pathlib import Path

# Extract chromosomes from appropriate file list
pattern = r'\.(chr[0-9X]+)\.'
if input_format == "plink":
    matches = [re.search(pattern, s) for s in plink_bim_files]

chrs = set([m.group(1) if m else None for m in matches])

print("Available chromosomes:", chrs)

# Load gene-trait pairs from CSV
gene_trait_pairs_df = pd.read_csv(gene_trait_pairs_to_test)

# Filter to the collection that are in the listed input files
gene_trait_pairs_df['chr'] = gene_trait_pairs_df['Region'].apply(get_gene_chr)
gene_trait_pairs_df = gene_trait_pairs_df[gene_trait_pairs_df['chr'].isin(chrs)]

# Extract unique genes and traits
genes = set(gene_trait_pairs_df['Region'])  # Second column is gene
traits = set(gene_trait_pairs_df['phenotype'])  # First column is trait

# Convert to sorted lists for consistency
genes = sorted(genes)
traits = sorted(traits)

print(f"Loaded {len(genes)} stable gene ID's.")
print(f"Stable gene IDs: {genes}")

# Load the phenotype data from JSON and store phenotype_IDs
with open(phenotype_json) as f:
    phenotype_data = json.load(f)

phenotype_ids = {entry["phenotype_ID"] for entry in phenotype_data}
print(f"Loaded {len(phenotype_ids)} phenotype_IDs from JSON.")
print(f"Phenotype IDs from JSON: {phenotype_ids}")

# Read model and variance ratio file lists
with open(list_of_loco_files) as f:
    loco_files = [line.strip() for line in f]

print("Loco Files Sample:", loco_files[:5])
print(f"All traits before filtering: {traits}")

def find_sep(infile):
    with open(infile) as f:
        headers = f.readline()
    if '\t' in headers:
        sep = '\t'
        engine=None
    elif ',' in headers:
        sep = ','
        engine=None
    else:
        sep = r'\s+'  # regex for any whitespace
        engine='python'
    return sep, engine

def set_fid_zero_inplace(infile):
    infile = Path(infile)
    suffix = infile.suffix.lower()
    
    # fam files: no header, always 6 columns
    if suffix in (".fam", ".id"):
        sep, engine = find_sep(infile)
        df = pd.read_csv(infile, sep=sep, engine=engine, header=None, dtype=str)
        df[0] = "0"
        df.to_csv(infile, sep=" ", header=False, index=False)
        print(f"✔ Overwrote {infile} with FID=0 for {len(df)} rows")
    
    # covariate and phenotype files: have headers
    else:
        sep, engine = find_sep(infile)
        df = pd.read_csv(infile, sep=sep, engine=engine, dtype='str')
        if "FID" not in df.columns or "IID" not in df.columns:
            raise ValueError(f"{infile} must contain 'FID' and 'IID' columns")
        df["FID"] = "0"
        df.to_csv(infile, sep="\t", header=True, index=False)
        print(f"✔ Overwrote {infile} with FID=0 for {len(df)} rows")

# Check and ensure that FID, IID in the phenotypes and covariates match the .fam file
# We assume that FID is set to 0 for all samples by default
set_fid_zero_inplace(covariate_file)
set_fid_zero_inplace(phenotype_file)
set_fid_zero_inplace(keep)

sep, engine = find_sep(phenotype_file)
headers = pd.read_csv(phenotype_file, sep=sep,
    engine=engine, dtype=str, nrows=0).columns.tolist()

# Filter traits based on presence in model and variance ratio files
available_traits = set()
for pid in phenotype_ids:  # phenotype IDs from JSON
    # Match start of filename or bounded by separators (_ . -)
    pattern = rf'(?:^|[/_.\-]){re.escape(pid)}(?=$|[/_.\-])'
    trait_in_loco = any(re.search(pattern, lf) for lf in loco_files)
    trait_in_phenotypes = any(re.search(pattern, pf) for pf in headers)
    if trait_in_loco and trait_in_phenotypes:
        available_traits.add(pid)

# For all available traits, determine the loco file and the phenotype, and create the
# pred file
run_files_dir = Path("run_files")
run_files_dir.mkdir(exist_ok=True)

for trait in available_traits:
    pattern = rf'(?:^|[/_.\-]){re.escape(trait)}(?=$|[/_.\-])'
    matching_files = [lf for lf in loco_files if re.search(pattern, lf)]
    if len(matching_files) != 1:
        raise ValueError(f"Trait '{trait}' expected exactly 1 match in loco_files, but found {len(matching_files)}: {matching_files}")
    loco_path = matching_files[0]
    # Write a separate TSV file for this trait
    output_file = f"run_files/{trait}_pred.list"
    with open(output_file, "w") as f:
        f.write(f"{trait} {loco_path}\n")

print(f"Wrote {len(available_traits)} space-delimited pred files to run_files")

# Store valid gene-trait pairs as a list
valid_gene_trait_pairs = [f"{gene}_{trait}" for gene, trait in zip(gene_trait_pairs_df.iloc[:, 1], gene_trait_pairs_df.iloc[:, 0]) if trait in available_traits]

print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

genes_in_valid_pairs = sorted({pair.split("_")[0] for pair in valid_gene_trait_pairs})

# Define format-specific output files for the 'all' rule
def get_format_outputs():
    base_outputs = [
        expand("run_files/{gene_trait}_{distance}_{maf}_string.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}.annotation.txt", gene=genes_in_valid_pairs),
        expand("run_files/{gene}.setlist.txt", gene=genes_in_valid_pairs),
        expand("run_files/bed/{gene}.bed", gene=genes_in_valid_pairs),
        expand("run_files/bed/expanded_regions_{gene}.bed", gene=genes_in_valid_pairs),
        expand("run_files/bed/expanded_coding_regions_{gene}.bed", gene=genes_in_valid_pairs),
        # expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               # gene_trait=valid_gene_trait_pairs,
               # distance=config["distance"],
               # maf=config["maf"]),
        # "brava_stepwise_conditional_analysis_results.txt"
    ]
    
    if input_format == "plink":
        base_outputs.extend([
            expand("run_files/{gene}_{distance}.bim", 
                   gene=genes_in_valid_pairs, distance=config["distance"]),
            expand("run_files/{gene}_{distance}.bed", 
                   gene=genes_in_valid_pairs, distance=config["distance"]),
            expand("run_files/{gene}_{distance}.fam", 
                   gene=genes_in_valid_pairs, distance=config["distance"]),
            expand("run_files/{gene}_{distance}_{maf}.bim", 
                   gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
            expand("run_files/{gene}_{distance}_{maf}.bed", 
                   gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
            expand("run_files/{gene}_{distance}_{maf}.fam", 
                   gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"])
        ])
    
    return base_outputs

# Target Rule for Completion of Pipeline
rule all:
    input:
        get_format_outputs()

rule identify_gene_start_stop:
    output:
        "run_files/bed/{gene}.bed",
        "run_files/bed/expanded_regions_{gene}.bed",
        "run_files/bed/expanded_coding_regions_{gene}.bed"
    params:
        distance=distance,
        outfolder="run_files/bed"
    log:
        stdout="logs/identify_gene_start_stop/{gene}.out",
        stderr="logs/identify_gene_start_stop/{gene}.err"
    shell:
        """
        set -euo pipefail
        python scripts/start_end_query.py \
            --ensembl_id \"{wildcards.gene}\" > {log.stdout} 2> {log.stderr}
        bash scripts/expand_coding_region.sh {wildcards.gene} \
            {params.distance} {params.outfolder} \
            > >(tee -a {log.stdout}) \
            2> >(tee -a {log.stderr} >&2)
        """

# PLINK-specific rules 
rule filter_to_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files if input_format == "plink" else [],
        plink_bed = lambda wildcards: plink_bed_files if input_format == "plink" else [],
        plink_fam = lambda wildcards: plink_fam_files if input_format == "plink" else [],
        keep = keep,
        regions = "run_files/bed/expanded_regions_{gene}.bed"
    output:
        bim = "run_files/{gene}_{distance}.bim",
        bed = "run_files/{gene}_{distance}.bed",
        fam = "run_files/{gene}_{distance}.fam"
    params:
        distance=distance,
        threads=config["threads"],
        outfolder="run_files"
    log:
        stdout="logs/filter_to_gene_plink/{gene}_{distance}.out",
        stderr="logs/filter_to_gene_plink/{gene}_{distance}.err"
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        echo $chr
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                matched_plink=$plink_fileset
                bash scripts/filter_to_gene_plink.sh $plink_fileset {wildcards.gene} \
                    {params.distance} {params.threads} {input.keep} \
                    {params.outfolder} \
                    > >(tee -a {log.stdout}) \
                    2> >(tee -a {log.stderr} >&2)
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
            exit 1
        fi
        """

rule filter_to_coding_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files if input_format == "plink" else [],
        plink_bed = lambda wildcards: plink_bed_files if input_format == "plink" else [],
        plink_fam = lambda wildcards: plink_fam_files if input_format == "plink" else [],
        keep = keep,
        regions = "run_files/bed/expanded_coding_regions_{gene}.bed"
    output:
        bim = "run_files/{gene}_{distance}_{maf}.bim",
        bed = "run_files/{gene}_{distance}_{maf}.bed",
        fam = "run_files/{gene}_{distance}_{maf}.fam"
    params:
        distance=distance,
        threads=config["threads"]
    log:
        stdout="logs/filter_to_coding_gene_plink/{gene}_{distance}_{maf}.out",
        stderr="logs/filter_to_coding_gene_plink/{gene}_{distance}_{maf}.err"
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        echo $chr
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                matched_plink=$plink_fileset
                bash scripts/filter_to_coding_gene_plink.sh $plink_fileset {wildcards.gene} \
                    {params.distance} {wildcards.maf} {params.threads} \
                    {input.keep} \
                    > >(tee -a {log.stdout}) \
                    2> >(tee -a {log.stderr} >&2)
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
            exit 1
        fi
        """

# Need to do this, looping over the files - figure out how to pass a list of files
rule filter_group_file_regenie:
    input:
        annotation = lambda wildcards: regenie_annotation_files,
        setlist = lambda wildcards: regenie_setlist_files
    output:
        annotation = "run_files/{gene}.annotation.txt",
        setlist = "run_files/{gene}.setlist.txt"
    log:
        stdout="logs/filter_group_file/{gene}.out",
        stderr="logs/filter_group_file/{gene}.err"
    params:
        out_prefix="run_files/{gene}",
        annotation_csv=lambda wildcards, input: ",".join(input.annotation),
        setlist_csv=lambda wildcards, input: ",".join(input.setlist)
    shell:
        """
        set -euo pipefail
        echo {params.annotation_csv}
        echo {params.setlist_csv}
        Rscript scripts/filter_group_file_regenie.R {wildcards.gene} {params.out_prefix} \
             {params.annotation_csv} {params.setlist_csv} \
            > >(tee -a {log.stdout}) \
            2> >(tee -a {log.stderr} >&2)
        """

# Format-specific stepwise conditional rules
rule spa_tests_stepwise_conditional_plink:
    input:
        pred_file = "run_files/{trait}_pred.list",
        phenotype_file = phenotype_file,
        covariate_file = covariate_file,
        plink_bim = "run_files/{gene}_{distance}_{maf}.bim",
        plink_bed = "run_files/{gene}_{distance}_{maf}.bed",
        plink_fam = "run_files/{gene}_{distance}_{maf}.fam",
        loco_file=lambda wildcards: [
            lf for lf in loco_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', lf)
        ],
        annotation_file = "run_files/{gene}.annotation.txt",
        setlist_file = "run_files/{gene}.setlist.txt"
    output:
        "run_files/{gene}_{trait}_{distance}_{maf}_string.txt"
    params:
        maf_common="{maf}",
        P_T=config["conditioning_pvalue"],
        covariate_cols = config["covariate_cols"],
        categ_covariate_cols = config["categorical_covariate_cols"]
    log:
        stdout="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.out",
        stderr="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.err"
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
            conda run --no-capture-output -n regenie_env \
                bash scripts/stepwise_conditional_regenie_plink.sh \
                $plink_fileset {output} {input.phenotype_file} {input.covariate_file} \
                {input.pred_file} $chr {params.P_T} \
                > >(tee -a {log.stdout}) \
                2> >(tee -a {log.stderr} >&2)
        done
        """

# # Format-specific conditional analysis rules
# rule spa_tests_conditional_plink:
#     input:
#         plink_bim = "run_files/{gene}_{distance}.bim",
#         plink_bed = "run_files/{gene}_{distance}.bed",
#         plink_fam = "run_files/{gene}_{distance}.fam",
#         model_file=lambda wildcards: [
#             mf for mf in model_files
#             if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', mf)
#         ],
#         variance_file=lambda wildcards: [
#             vf for vf in variance_files
#             if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', vf)
#         ],
#         sparse_matrix=sparse_matrix,
#         group_file="run_files/{gene}_group_file.txt",
#         conditioning_variants="run_files/{gene}_{trait}_{distance}_{maf}_string.txt"
#     output:
#         "saige_outputs/{gene}_{trait}_{distance}_saige_results_{maf}.txt" 
#     params:
#         min_mac=min_mac,
#         annotations_to_include=annotations_to_include,
#         max_MAF="{maf}",
#         use_null_var_ratio=config["use_null_var_ratio"]
#     log:
#         stdout="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.out",
#         stderr="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.err"
#     threads: 1
#     shell:
#         """
#         set -euo pipefail
#         plink_fileset=$(echo {input.plink_bed} | sed 's/\\.bed$//')
#         conda run --no-capture-output -n RSAIGE_vcf_version \
#             bash scripts/saige_step2_conditioning_check_plink.sh \
#             $plink_fileset {output} {params.min_mac} {input.model_file} \
#             {input.variance_file} {input.sparse_matrix} {input.group_file} \
#             {params.annotations_to_include} {input.conditioning_variants} \
#             {params.max_MAF} {params.use_null_var_ratio} \
#             > >(tee -a {log.stdout}) \
#             2> >(tee -a {log.stderr} >&2)
#         """

# rule combine_results:
#     input:
#         expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
#                gene_trait=valid_gene_trait_pairs,
#                distance=config["distance"],
#                maf=config["maf"]),
#     output:
#         "brava_stepwise_conditional_analysis_results.txt"
#     log:
#         stdout="logs/combine_results/final_output.out",
#         stderr="logs/combine_results/final_output.err"
#     shell:
#         """
#         set -euo pipefail
#         python scripts/combine_saige_outputs.py --out {output} \
#         > >(tee -a {log.stdout}) \
#         2> >(tee -a {log.stderr} >&2)
#         """

# # Rule order to ensure proper execution based on input format
ruleorder: filter_group_file_regenie > identify_gene_start_stop
# ruleorder: identify_gene_start_stop > filter_to_coding_gene_plink
# ruleorder: spa_tests_stepwise_conditional_plink > filter_to_gene_plink
