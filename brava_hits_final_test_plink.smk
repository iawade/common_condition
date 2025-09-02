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

min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re
from scripts.extract_chromosome import get_gene_chr

pattern = r'\.(chr[0-9X]+)\.'
matches = [re.search(pattern, s) for s in plink_bim_files]
chrs = set([m.group(1) if m else None for m in matches])

print("Available chromosomes:", chrs)

df = pd.read_csv(gene_trait_pairs_to_test_with_conditioning_variants, sep='\t')
df['chr'] = df['Gene'].apply(get_gene_chr)
df = df[df['chr'].isin(chrs)]

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
    filename = f"final_run_files/{job['Gene']}_{job['Trait']}_{job['MAF_cutoff_for_conditioning_variants']}_extract.bed"
    with open(filename, "w") as f:
        variants = job['cond'].split(",")
        for v in variants:
            chrom, pos, ref, alt = v.split(":")
            _ = f.write(f"{chrom}\t{pos}\t{pos}\n")

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
        "final_run_files/{gene}_{trait}_{maf}_ld_pruned_string.txt" 
    params:
        file="final_run_files/{gene}_{trait}_{maf}_ld_pruned_string"
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            if [[ "$plink_bed" =~ \\.($chr)\\. ]]; then
                plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
                matched_plink=$plink_fileset
                plink2 --bfile $plink_fileset \
                  --extract {input.conditioning_variants} \
                  --indep-pairwise 50 5 0.9 \
                  --out {params.file} || true
                if [[ -f {params.file}.prune.in ]]; then
                    paste -sd, {params.file}.prune.in > {output}
                else
                    # No variants to prune, create an empty output
                    touch {output}
                fi
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
        fi
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
        conditioning_variants="final_run_files/{gene}_{trait}_{distance}_{maf}_ld_pruned_string.txt"
    output:
        "final_saige_outputs/{gene}_{trait}_{distance}_saige_conditioned_results_{maf}.txt" 
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    threads: 8
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
        expand("final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}.txt",
        zip,
        gene=[job['Gene'] for job in conditioning_jobs],
        trait=[job['Trait'] for job in conditioning_jobs],
        maf=[job['MAF_cutoff_for_conditioning_variants'] for job in conditioning_jobs]
        ),
    output:
        "brava_final_conditional_analysis_results.txt"
    shell:
        """
        set -euo pipefail
        python scripts/combine_saige_outputs.py --out {output}
        """
