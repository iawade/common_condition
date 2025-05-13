# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis

configfile: "config.yaml"

# Read inputs from config file
list_of_vcf_files = config["list_of_vcf_files"]
list_of_model_files = config["list_of_model_files"]
list_of_variance_ratio_files = config["list_of_variance_ratio_files"]
list_of_group_files = config["list_of_group_files"]
sparse_matrix = config["sparse_matrix"]
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
phenotype_json = config["phenotype_json"]

# Load VCF files
with open(list_of_vcf_files) as f:
    vcf_files = [line.strip() for line in f]

print("VCF files:", vcf_files[:5])

# Load group files
with open(list_of_group_files) as f:
    group_files = [line.strip() for line in f]

print("Group files:", group_files[:5])

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

print(f"Loaded {len(genes)} stable gene ID's.")

# Debugging: Print phenotype IDs
print(f"Stable gene IDs from JSON: {genes}")

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

with open(list_of_variance_ratio_files) as f:
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

# Store valid gene-trait pairs as a list
valid_gene_trait_pairs = [f"{gene}_{trait}" for gene, trait in zip(gene_trait_pairs_df.iloc[:, 0], gene_trait_pairs_df.iloc[:, 1]) if trait in available_traits]

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

# Target Rule for Completion of Pipeline
rule all:
    input:
        # expand("run_files/{gene}_{distance}_{maf}_string.txt",
        # gene=genes, distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_group_file.txt", gene=genes) #,
        # expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               # gene_trait=valid_gene_trait_pairs,
               # distance=config["distance"],
               # maf=config["maf"]),
        # "brava_conditional_analysis_results.txt"

rule identify_gene_start_stop:
    output:
        "run_files/{gene}.bed"
    shell:
        "python scripts/start_end_query.py --ensembl_id \"{wildcards.gene}\""

rule id_variants_for_conditioning:
    input:
        vcf = lambda wildcards: vcf_files,
        bed = "run_files/{gene}.bed" 
    output:
        "run_files/{gene}_{distance}_{maf}_list.txt",
        "run_files/{gene}_{distance}_{maf}_string.txt"
    params:
        distance=distance,
        threads=config["threads"]
    shell:
        """
        for vcf in {input.vcf}; do
            bash scripts/id_variants_for_common_variant_conditioning.sh $vcf {wildcards.gene} {params.distance} {wildcards.maf} {params.threads}
        done
        """

rule filter_group_file:
    input:
        group = lambda wildcards: group_files
    output:
        "run_files/{gene}_group_file.txt"
    shell:
        """
        for group in {input.group}; do
            if [[ "$group" == *.gz ]]; then
                   zcat "$group" | grep -m1 -A1 "{wildcards.gene}" >> {output}
            else
                   grep -m1 -A1 "{wildcards.gene}" "$group" >> {output}
            fi
        done
        touch {output}
        """

rule spa_tests_conditional:
    input:
        vcf=lambda wildcards: vcf_files,
        model_file=lambda wildcards: [mf for mf in model_files if wildcards.trait in mf],  
        variance_file=lambda wildcards: [vf for vf in variance_files if wildcards.trait in vf],    
        sparse_matrix=sparse_matrix,
        group_file="run_files/{gene}_group_file.txt",
        conditioning_variants="run_files/{gene}_{distance}_{maf}_string.txt"
    output:
        "saige_outputs/{gene}_{trait}_{distance}_saige_results_{maf}.txt" 
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}"
    shell:
        """
        for vcf in {input.vcf}; do
            bash scripts/saige_step2_conditioning_check.sh \
                $vcf {output} {params.min_mac} {input.model_file} {input.variance_file} {input.sparse_matrix} {input.group_file} {params.annotations_to_include} {input.conditioning_variants} {params.max_MAF}
        done
        """

rule combine_results:
    input:
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
    output:
        "brava_conditional_analysis_results.txt",
    shell:
        """
        python scripts/combine_saige_outputs.py --out {output}
        """
