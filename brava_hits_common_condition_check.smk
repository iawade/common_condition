# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditioninal Analysis

configfile: "config.yaml"

# Read inputs from config file
list_of_vcf_files = config["list_of_vcf_files"]
list_of_model_files = config["list_of_model_files"]
variance_ratio_file = config["variance_ratio_file"]
sparse_matrix = config["sparse_matrix"]
group_file = config["group_file"]
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
phenotype_json = config["phenotype_json"]

# Load VCF files
with open(list_of_vcf_files) as f:
    vcf_files = [line.strip() for line in f]

# Read model and variance ratio files
with open(list_of_model_files) as f:
    model_files = [line.strip() for line in f]

with open(variance_ratio_file) as f:
    variance_files = [line.strip() for line in f]

distance = config["distance"]
maf = config["maf"]

min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json

# Load gene-trait pairs from CSV
gene_trait_pairs_df = pd.read_csv(gene_trait_pairs_to_test)

# Extract unique genes and traits
genes = set(gene_trait_pairs_df.iloc[:, 0])  # First column as gene
traits = set(gene_trait_pairs_df.iloc[:, 1])  # Second column as trait

# Convert to sorted lists for consistency
genes = sorted(genes)
traits = sorted(traits)

# Load the phenotype data from JSON and store phenotype_IDs
with open(phenotype_json) as f:
    phenotype_data = json.load(f)

phenotype_ids = {entry["phenotype_ID"] for entry in phenotype_data}
print(f"Loaded {len(phenotype_ids)} phenotype_IDs from JSON.")

# Debugging: Print phenotype IDs
print(f"Phenotype IDs from JSON: {phenotype_ids}")

# Read model and variance ratio file lists
with open(list_of_model_files) as f:
    model_files = [line.strip() for line in f]

with open(variance_ratio_file) as f:
    variance_files = [line.strip() for line in f]

# Debugging: Print first few model and variance files
print("Model Files Sample:", model_files[:5])
print("Variance Files Sample:", variance_files[:5])

# Debugging: Print traits before filtering
print(f"All traits before filtering: {traits}")

# Filter traits based on presence in model and variance ratio files
available_traits = set()
for pid in phenotype_ids:  # phenotype IDs from JSON
    trait_in_model = any(pid in mf for mf in model_files)
    trait_in_variance = any(pid in vf for vf in variance_files)
    
    if trait_in_model and trait_in_variance:
        available_traits.add(pid)  # Store the phenotype ID instead of an incorrect trait name

# Debugging: Print traits that passed filtering
print(f"Available traits after filtering: {available_traits}")

# Store valid gene-trait pairs as tuples
valid_gene_trait_pairs = {(gene, trait) for gene, trait in zip(gene_trait_pairs_df.iloc[:, 0], gene_trait_pairs_df.iloc[:, 1]) if trait in available_traits}

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

# Debugging: Print traits that passed filtering
print(f"Available traits after filtering: {available_traits}")

# Store valid gene-trait pairs as tuples
valid_gene_trait_pairs = {(gene, trait) for gene, trait in zip(gene_trait_pairs_df.iloc[:, 0], gene_trait_pairs_df.iloc[:, 1]) if trait in available_traits}

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

print(f"Loaded {len(valid_gene_trait_pairs)} gene-trait pairs.")

# Define chromosome list
chromosomes = [f"chr{i}" for i in range(1, 23)]
print(f"Chromosomes: {chromosomes}")

# Target Rule for Completion of Pipeline
rule all:
    input:
        "brava_conditional_analysis_results.txt"

rule identify_gene_start_stop:
    output:
        "{gene}.bed"
    shell:
        "python scripts/biomart_start_end_query.py --ensembl_id \"{wildcards.gene}\""

rule id_variants_for_conditioning:
    input:
        lambda wildcards: vcf_files,
        "{gene}.bed"
    output:
        "{gene}_{distance}_{maf}_list.txt",
        "{gene}_{distance}_{maf}_string.txt"
    params:
        distance=distance,
        threads=config["threads"]
    shell:
        """
        bash scripts/id_variants_for_common_variant_conditioning.sh {input[0]} {wildcards.gene} {params.distance} {wildcards.maf} {params.threads}
        """

rule filter_group_file:
    input:
        group_file,
    output:
        "{gene}_group_file.txt",
    shell:
        """
        grep {wildcards.gene} {input[0]} > {output} || touch {output}
        """

rule spa_tests_conditional:
    input:
        lambda wildcards: vcf_files,
        lambda wildcards: [mf for mf in model_files if wildcards.trait in mf],
        lambda wildcards: [vf for vf in variance_files if wildcards.trait in vf],
        sparse_matrix,
        "{gene}_group_file.txt",
        lambda wildcards: f"{wildcards.gene}_{distance}_{wildcards.maf}_string.txt"
    output:
        "{gene}_{trait}_{chrom}_saige_results_{maf}.txt"
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include
    shell:
        """
        bash scripts/saige_step2_conditioning_check.sh \
         {input[0]} {output} {wildcards.chrom} {params.min_mac} {input[1]} {input[2]} {input[3]} {input[4]} {params.annotations_to_include} {input[5]}
        """

rule combine_results:
    input:
        expand(
            "{gene}_{trait}_{chrom}_saige_results_{maf}.txt",
            gene=[gene for gene, trait in valid_gene_trait_pairs],  # Use only valid genes
            trait=[trait for gene, trait in valid_gene_trait_pairs],  # Use only valid traits
            chrom=chromosomes,
            maf=config["maf"]
        )
    output:
        "brava_conditional_analysis_results.txt",
    params:
    shell:
        """
        cat {input} > {output}
        """