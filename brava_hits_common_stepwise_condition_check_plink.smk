# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis

configfile: "config.yaml"

# Read inputs from config file
list_of_plink_files = config["list_of_plink_files"]
list_of_model_files = config["list_of_model_files"]
list_of_variance_ratio_files = config["list_of_variance_ratio_files"]
list_of_group_files = config["list_of_group_files"]
sparse_matrix = config["sparse_matrix"]
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
phenotype_json = config["phenotype_json"]
use_null_var_ratio = config["use_null_var_ratio"]

# For debugging
# list_of_plink_files="plink_testing/afr_plink_file_list.txt"
# list_of_model_files="plink_testing/afr_model_file_list.txt"
# list_of_variance_ratio_files="plink_testing/afr_variance_file_list.txt"
# sparse_matrix="/home/jupyter/conditioning_test/afr/mtx/allofus_array_afr_snp_wise_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
# list_of_group_files="plink_testing/afr_group_file_list.txt"
# distance=500000
# maf=[0.001, 0.0001]
# min_mac=0.5
# annotations_to_include="pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous"
# gene_trait_pairs_to_test="gene_pheno_test.csv"
# phenotype_json="pilot_phenotypes.json"

# Load plink files
with open(list_of_plink_files) as f:
    plink_files = [line.strip() for line in f]

plink_bim_files = [f"{p}.bim" for p in plink_files]
plink_bed_files = [f"{p}.bed" for p in plink_files]
plink_fam_files = [f"{p}.fam" for p in plink_files]

# Load group files
with open(list_of_group_files) as f:
    group_files = [line.strip() for line in f]

print("Plink Files Sample:", plink_files[:5])
print("Group Files Sample:", group_files[:5])

distance = config["distance"]
maf = config["maf"]

min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re

# Load gene-trait pairs from CSV
gene_trait_pairs_df = pd.read_csv(gene_trait_pairs_to_test)

# Extract unique genes and traits
genes = set(gene_trait_pairs_df.iloc[:, 1])  # Second column is gene
traits = set(gene_trait_pairs_df.iloc[:, 0])  # First column is trait

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
    # Match start of filename or bounded by separators (_ . -)
    pattern = rf'(?:^|[/_.\-]){re.escape(pid)}(?=[/_.\-])'
    trait_in_model = any(re.search(pattern, mf) for mf in model_files)
    trait_in_variance = any(re.search(pattern, mf) for mf in model_files)
    if trait_in_model and trait_in_variance:
        available_traits.add(pid)

# Store valid gene-trait pairs as a list
valid_gene_trait_pairs = [f"{gene}_{trait}" for gene, trait in zip(gene_trait_pairs_df.iloc[:, 1], gene_trait_pairs_df.iloc[:, 0]) if trait in available_traits]

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

genes_in_valid_pairs = sorted({pair.split("_")[0] for pair in valid_gene_trait_pairs})

# Target Rule for Completion of Pipeline
rule all:
    input:
        expand("run_files/{gene_trait}_{distance}_{maf}_string.txt",
        gene_trait=valid_gene_trait_pairs,
        distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_{distance}_{maf}.bim", 
        gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_{distance}_{maf}.bed", 
        gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_{distance}_{maf}.fam", 
        gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_group_file.txt", gene=genes_in_valid_pairs),
        expand("run_files/{gene}.bed", gene=genes_in_valid_pairs),
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
        "brava_stepwise_conditional_analysis_results.txt"

rule identify_gene_start_stop:
    output:
        "run_files/start_end_{gene}.bed"
    shell:
        "python scripts/start_end_query.py --ensembl_id \"{wildcards.gene}\""

rule filter_to_coding_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files,
        plink_bed = lambda wildcards: plink_bed_files,
        plink_fam = lambda wildcards: plink_fam_files,
        bed = "run_files/{gene}.bed" 
    output:
        "run_files/{gene}_{distance}_{maf}.bim",
        "run_files/{gene}_{distance}_{maf}.bed",
        "run_files/{gene}_{distance}_{maf}.fam"
    params:
        distance=distance,
        threads=config["threads"]
    shell:
        """
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        echo $chr
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                matched_plink=$plink_fileset
                bash scripts/filter_to_coding_gene_plink.sh $plink_fileset {wildcards.gene} {params.distance} {wildcards.maf} {params.threads}
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
        fi
        """

rule filter_group_file:
    input:
        group = lambda wildcards: group_files
    output:
        "run_files/{gene}_group_file.txt"
    shell:
        """
        > {output}
        for group in {input.group}; do
            if [[ "$group" == *.gz ]]; then
                   zcat "$group" | grep -m1 -A1 "{wildcards.gene}" >> {output} || true
            else
                   grep -m1 -A1 "{wildcards.gene}" "$group" >> {output} || true
            fi
        done
        touch {output}
        """

rule spa_tests_stepwise_conditional:
    input:
        plink_bim = "run_files/{gene}_{distance}_{maf}.bim",
        plink_bed = "run_files/{gene}_{distance}_{maf}.bed",
        plink_fam = "run_files/{gene}_{distance}_{maf}.fam",
        model_file=lambda wildcards: [
            mf for mf in model_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', mf)
        ],
        variance_file=lambda wildcards: [
            vf for vf in variance_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', vf)
        ],
        sparse_matrix=sparse_matrix,
        group_file="run_files/{gene}_group_file.txt",
    output:
        "run_files/{gene}_{trait}_{distance}_{maf}_string.txt" 
    params:
        maf_common="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    shell:
        """
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
            bash scripts/stepwise_conditional_SAIGE_plink.sh \
                $plink_fileset {output} {input.model_file} {input.variance_file} {input.sparse_matrix} $chr {params.use_null_var_ratio}
        done
        """

rule spa_tests_conditional:
    input:
        plink_bim = lambda wildcards: plink_bim_files,
        plink_bed = lambda wildcards: plink_bed_files,
        plink_fam = lambda wildcards: plink_fam_files,
        model_file=lambda wildcards: [
            mf for mf in model_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', mf)
        ],
        variance_file=lambda wildcards: [
            vf for vf in variance_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', vf)
        ],
        sparse_matrix=sparse_matrix,
        group_file="run_files/{gene}_group_file.txt",
        conditioning_variants="run_files/{gene}_{trait}_{distance}_{maf}_string.txt"
    output:
        "saige_outputs/{gene}_{trait}_{distance}_saige_results_{maf}.txt" 
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    shell:
        """
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                bash scripts/saige_step2_conditioning_check_plink.sh \
                    $plink_fileset {output} {params.min_mac} {input.model_file} {input.variance_file} {input.sparse_matrix} {input.group_file} {params.annotations_to_include} {input.conditioning_variants} {params.max_MAF} {params.use_null_var_ratio}
            fi
        done
        """

rule combine_results:
    input:
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
    output:
        "brava_stepwise_conditional_analysis_results.txt"
    shell:
        """
        python scripts/combine_saige_outputs.py --out {output}
        """
