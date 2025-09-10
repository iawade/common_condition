# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis

configfile: "config.yaml"

# Read inputs from config file
list_of_plink_files = config["list_of_plink_files"]
list_of_model_files = config["list_of_model_files"]
list_of_variance_ratio_files = config["list_of_variance_ratio_files"]
list_of_group_files = config["list_of_group_files"]
sparse_matrix = config["sparse_matrix"]
sparse_matrix_id = f"{sparse_matrix}.sampleIDs.txt"
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
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

print("Plink Files Sample:", plink_files[:5])
print("Group Files Sample:", group_files[:5])

distance = config["distance"]
maf = config["maf"]

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

# Load gene-trait pairs from CSV
gene_trait_pairs_df = pd.read_csv(gene_trait_pairs_to_test)

# Filter to the collection that are in the listed vcfs
gene_trait_pairs_df['chr'] = gene_trait_pairs_df['Region'].apply(get_gene_chr)
gene_trait_pairs_df = gene_trait_pairs_df[gene_trait_pairs_df['chr'].isin(chrs)]

# Extract unique genes and traits
genes = set(gene_trait_pairs_df['Region'])  # Second column is gene
traits = set(gene_trait_pairs_df['phenotype']) # First column is trait

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
        expand("run_files/bed/{gene}.bed", gene=genes_in_valid_pairs),
        expand("run_files/bed/expanded_regions_{gene}.bed", gene=genes_in_valid_pairs),
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
        "brava_stepwise_conditional_analysis_results.txt"

rule identify_gene_start_stop:
    output:
        r"run_files/bed/{gene,[^/]+}.bed",
        r"run_files/bed/expanded_regions_{gene,[^/]+}.bed"
    params:
        distance=distance
    shell:
        """
        set -euo pipefail
        python scripts/start_end_query.py --ensembl_id \"{wildcards.gene}\"
        bash scripts/expand_coding_region.sh {wildcards.gene} {params.distance}
        """

rule filter_to_coding_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files,
        plink_bed = lambda wildcards: plink_bed_files,
        plink_fam = lambda wildcards: plink_fam_files,
        sparse_matrix_id = sparse_matrix_id,
        regions = r"run_files/bed/expanded_regions_{gene,[^/]+}.bed"
    output:
        bim = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.bim",
        bed = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.bed",
        fam = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.fam"
    params:
        distance=distance,
        threads=config["threads"]
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
                bash scripts/filter_to_coding_gene_plink.sh $plink_fileset {wildcards.gene} {params.distance} {wildcards.maf} {params.threads} {input.sparse_matrix_id}
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

rule spa_tests_stepwise_conditional:
    input:
        plink_bim = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.bim",
        plink_bed = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.bed",
        plink_fam = r"run_files/{gene,[^/]+}_{distance,\d+}_{maf,[0-9.]+}.fam",
        model_file=lambda wildcards: [
            mf for mf in model_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', mf)
        ],
        variance_file=lambda wildcards: [
            vf for vf in variance_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', vf)
        ],
        sparse_matrix=sparse_matrix,
        sparse_matrix_id=sparse_matrix_id,
        group_file="run_files/{gene}_group_file.txt",
    output:
        "run_files/{gene}_{trait}_{distance}_{maf}_string.txt" 
    params:
        maf_common="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
            conda run --no-capture-output --prefix envs/RSAIGE_vcf_version bash scripts/stepwise_conditional_SAIGE_plink.sh \
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
        "brava_stepwise_conditional_analysis_results.txt"
    shell:
        """
        set -euo pipefail
        python scripts/combine_saige_outputs.py --out {output}
        """

ruleorder:
    identify_gene_start_stop > filter_to_coding_gene_plink
