# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis

configfile: "config.yaml"

# Read inputs from config file
list_of_plink_files = config["list_of_plink_files"]
list_of_model_files = config["list_of_model_files"]
list_of_variance_ratio_files = config["list_of_variance_ratio_files"]
list_of_group_files = config["list_of_group_files"]
sparse_matrix = config["sparse_matrix"]
gene_trait_pairs_to_test_with_conditioning_variants = config["gene_trait_pairs_to_test_with_conditioning_variants"]
phenotype_json = config["phenotype_json"]
use_null_var_ratio = config["use_null_var_ratio"]

# Load plink files
with open(list_of_plink_files) as f:
    plink_files = [line.strip() for line in f]

plink_bim_files = [f"{p}.bim" for p in plink_files]
plink_bed_files = [f"{p}.bed" for p in plink_files]
plink_fam_files = [f"{p}.fam" for p in plink_files]

# Load group files
with open(list_of_group_files) as f:
    group_files = [line.strip() for line in f]

distance = config["distance"]
maf = config["maf"]

min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re

df = pd.read_csv("/Users/dpalmer/Repositories/BRaVa_curation/meta_analysis/conditional_results/data/gene_phenotype_pairs_with_conditioning_variants.tsv", sep='\t')
conditioning_jobs = df.to_dict(orient='records')

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

conditioning_jobs = [job for job in conditioning_jobs if job['Trait'] in available_traits]

for job in conditioning_jobs:
    filename = f"final_run_files/{job['Gene']}_{job['Trait']}_{job['MAF_cutoff_for_conditioning_variants']}_extract.txt"
    with open(filename, "w") as f:
        f.writelines(f"{v}\n" for v in job['cond'].split(","))

genes = sorted(set(job['Gene'] for job in conditioning_jobs))

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(conditioning_jobs)} conditioning jobs with available model/variance files.")

# Target Rule for Completion of Pipeline
rule all:
    input:
        expand([
            "final_run_files/{gene}_{trait}_{maf}_extract.txt"
            "final_run_files/{gene}_{trait}_{maf}_ld_pruned_string.txt",
            "final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}.txt"],
        zip,
        gene=[job['Gene'] for job in conditioning_jobs],
        trait=[job['Trait'] for job in conditioning_jobs],
        maf=[job['MAF_cutoff_for_conditioning_variants'] for job in conditioning_jobs]
        ),
        expand("final_run_files/{gene}_group_file.txt", gene=genes),
        "brava_final_conditional_analysis_results.txt"

rule filter_group_file:
    input:
        group = lambda wildcards: group_files
    output:
        "final_run_files/{gene}_group_file.txt"
    shell:
        """
        set -euo pipefail
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

rule prune_to_independent_conditioning_variants:
    input:
        plink_bim = lambda wildcards: plink_bim_files,
        plink_bed = lambda wildcards: plink_bed_files,
        plink_fam = lambda wildcards: plink_fam_files,
        conditioning_variants = "final_run_files/{gene}_{trait}_{maf}_extract.txt",
        group_file="final_run_files/{gene}_group_file.txt",
    output:
        "final_run_files/{gene}_{trait}_{maf}_ld_pruned_string.prune.in" 
    params:
        file="final_run_files/{gene}_{trait}_{maf}_ld_pruned_string"
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                plink2 --bfile $plink_fileset \
                  --extract {input.conditioning_variants} \
                  --indep-pairwise 50 5 0.9 \
                  --out {params.file}
            fi
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
        group_file="final_run_files/{gene}_group_file.txt",
        conditioning_variants="final_run_files/{gene}_{trait}_{distance}_{maf}_string.txt"
    output:
        "final_saige_outputs/{gene}_{trait}_{distance}_saige_results_{maf}.txt" 
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    threads: 4
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                conda run --no-capture-output --prefix envs/RSAIGE_vcf_version bash scripts/saige_step2_conditioning_check_plink.sh \
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
        "brava_final_conditional_analysis_results.txt"
    shell:
        """
        set -euo pipefail
        python scripts/combine_saige_outputs.py --out {output}
        """
