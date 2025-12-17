# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis
# Supports both VCF and PLINK input formats

configfile: "config.yaml"

wildcard_constraints:
    gene="[A-Za-z0-9]+",
    trait="[A-Za-z0-9]+",
    gene_trait="[A-Za-z0-9]+_[A-Za-z0-9]+",  # For concatenated gene_trait wildcards
    distance="[0-9]+",
    maf="[0-9.]+"

# Read inputs from config file
input_format = config["input_format"]  # "vcf" or "plink"
sparse_matrix = config["sparse_matrix"]
sparse_matrix_id = f"{sparse_matrix}.sampleIDs.txt"
gene_trait_pairs_to_test = config["gene_trait_pairs_to_test"]
protein_coding_region_bed = config["protein_coding_region_bed"]
phenotype_json = config["phenotype_json"]
use_null_var_ratio = config["use_null_var_ratio"]

# Format-specific file lists
if input_format == "vcf":
    list_of_input_files = config["list_of_vcf_files"]
elif input_format == "plink":
    list_of_input_files = config["list_of_plink_files"]
else:
    raise ValueError(f"Invalid input_format: {input_format}. Must be 'vcf' or 'plink'")

list_of_model_files = config["list_of_model_files"]
list_of_variance_ratio_files = config["list_of_variance_ratio_files"]
list_of_group_files = config["list_of_group_files"]

# Load input files
with open(list_of_input_files) as f:
    input_files = [line.strip() for line in f]

# Create format-specific file lists for PLINK
if input_format == "plink":
    plink_bim_files = [f"{p}.bim" for p in input_files]
    plink_bed_files = [f"{p}.bed" for p in input_files]
    plink_fam_files = [f"{p}.fam" for p in input_files]

# Load group files
with open(list_of_group_files) as f:
    group_files = [line.strip() for line in f]

print(f"Input format: {input_format}")
print(f"Gene trait pairs to test: {gene_trait_pairs_to_test}")
print(f"Input files Sample: {input_files[:5]}")
print(f"Group files Sample: {group_files[:5]}")

distance = config["distance"]
maf = config["maf"]
min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re
from scripts.extract_chromosome import get_gene_chr

# Extract chromosomes from appropriate file list
pattern = r'\.(chr[0-9X]+)\.'
if input_format == "plink":
    matches = [re.search(pattern, s) for s in plink_bim_files]
else:  # vcf
    matches = [re.search(pattern, s) for s in input_files]
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
with open(list_of_model_files) as f:
    model_files = [line.strip() for line in f]

with open(list_of_variance_ratio_files) as f:
    variance_files = [line.strip() for line in f]

print("Model Files Sample:", model_files[:5])
print("Variance Files Sample:", variance_files[:5])
print(f"All traits before filtering: {traits}")

# Filter traits based on presence in model and variance ratio files
available_traits = set()
for pid in phenotype_ids:  # phenotype IDs from JSON
    # Match start of filename or bounded by separators (_ . -)
    pattern = rf'(?:^|[/_.\-]){re.escape(pid)}(?=[/_.\-])'
    trait_in_model = any(re.search(pattern, mf) for mf in model_files)
    trait_in_variance = any(re.search(pattern, vf) for vf in variance_files)
    if trait_in_model and trait_in_variance:
        available_traits.add(pid)

# Store valid gene-trait pairs as a list
valid_gene_trait_pairs = [f"{gene}_{trait}" for gene, trait in zip(gene_trait_pairs_df.iloc[:, 1], gene_trait_pairs_df.iloc[:, 0]) if trait in available_traits]

print(f"Filtered {len(valid_gene_trait_pairs)} gene-trait pairs with available model/variance files.")
print(f"Valid gene-trait pairs: {valid_gene_trait_pairs}")

genes_in_valid_pairs = sorted({pair.split("_")[0] for pair in valid_gene_trait_pairs})

# Group-file sanity check
present_genes = set()

for gf in group_files:
    try:
        with open(gf) as f:
            for line in f:
                if not line.strip():
                    continue
                # space-delimited: first column only
                gene_id = line.split()[0]
                present_genes.add(gene_id)
    except OSError:
        print(f"WARNING: Could not read group file {gf}")

genes_expected = set(genes_in_valid_pairs)
n_expected = len(genes_expected)
n_present = len(genes_expected & present_genes)
pct_present = (n_present / n_expected * 100) if n_expected else 0.0

print(
    f"Group-file presence check: "
    f"{n_present}/{n_expected} genes ({pct_present:.1f}%) present"
)

if pct_present < 90.0:
    missing = sorted(genes_expected - present_genes)
    print(
        "WARNING: More than 10% of genes are missing from group files "
        f"({len(missing)} missing)."
    )
    print("Example missing genes (up to 10):", missing[:10])

# Define format-specific output files for the 'all' rule
def get_format_outputs():
    base_outputs = [
        expand("run_files/{gene_trait}_{distance}_{maf}_string.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"], maf=config["maf"]),
        expand("run_files/{gene}_group_file.txt", gene=genes_in_valid_pairs),
        expand("run_files/bed/{gene}.bed", gene=genes_in_valid_pairs),
        expand("run_files/bed/expanded_regions_{gene}.bed", gene=genes_in_valid_pairs),
        expand("run_files/bed/expanded_coding_regions_{gene}.bed", gene=genes_in_valid_pairs),
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
        "brava_stepwise_conditional_analysis_results.txt"
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
    else:  # vcf
        base_outputs.extend([
            expand("run_files/{gene}_{distance}.vcf.bgz", 
                   gene=genes_in_valid_pairs, distance=config["distance"]),
            expand("run_files/{gene}_{distance}.vcf.bgz.csi", 
                   gene=genes_in_valid_pairs, distance=config["distance"]),
            expand("run_files/{gene}_{distance}_{maf}.vcf.bgz", 
                   gene=genes_in_valid_pairs, distance=config["distance"], maf=config["maf"]),
            expand("run_files/{gene}_{distance}_{maf}.vcf.bgz.csi", 
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

# VCF-specific rules
rule filter_to_gene_vcf:
    input:
        vcf = lambda wildcards: input_files if input_format == "vcf" else [],
        bed = "run_files/bed/expanded_regions_{gene}.bed"
    output:
        "run_files/{gene}_{distance}.vcf.bgz",
        "run_files/{gene}_{distance}.vcf.bgz.csi"
    params:
        distance=distance,
        threads=config["threads"],
        outfolder="run_files"
    log:
        stdout="logs/filter_to_gene_vcf/{gene}_{distance}.out",
        stderr="logs/filter_to_gene_vcf/{gene}_{distance}.err"
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        echo $chr
        for vcf in {input.vcf}; do
            if [[ "$vcf" =~ \\.($chr)\\. ]]; then
                matched_vcf=$vcf
                bash scripts/filter_to_gene_vcf.sh $vcf {wildcards.gene} \
                    {params.distance} {params.threads} \
                    {params.outfolder} \
                    > >(tee -a {log.stdout}) \
                    2> >(tee -a {log.stderr} >&2)
            fi
        done

        if [[ -z "$matched_vcf" ]]; then
            echo "ERROR: No matching VCF found for chromosome $chr"
            exit 1
        fi
        """

rule filter_to_coding_gene_vcf:
    input:
        vcf = lambda wildcards: input_files if input_format == "vcf" else [],
        bed = "run_files/bed/expanded_coding_regions_{gene}.bed" 
    output:
        "run_files/{gene}_{distance}_{maf}.vcf.bgz",
        "run_files/{gene}_{distance}_{maf}.vcf.bgz.csi"
    params:
        distance=distance,
        threads=config["threads"],
        outfolder="run_files"
    log:
        stdout="logs/filter_to_coding_gene_vcf/{gene}_{distance}_{maf}.out",
        stderr="logs/filter_to_coding_gene_vcf/{gene}_{distance}_{maf}.err"
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        echo $chr
        for vcf in {input.vcf}; do
            if [[ "$vcf" =~ \\.($chr)\\. ]]; then
                matched_vcf=$vcf
                bash scripts/filter_to_coding_gene_vcf.sh $vcf {wildcards.gene} \
                    {params.distance} {wildcards.maf} \
                    {params.threads} \
                    > >(tee -a {log.stdout}) \
                    2> >(tee -a {log.stderr} >&2)
            fi
        done

        if [[ -z "$matched_vcf" ]]; then
            echo "ERROR: No matching VCF found for chromosome $chr"
            exit 1
        fi
        """

# PLINK-specific rules 
rule filter_to_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files if input_format == "plink" else [],
        plink_bed = lambda wildcards: plink_bed_files if input_format == "plink" else [],
        plink_fam = lambda wildcards: plink_fam_files if input_format == "plink" else [],
        sparse_matrix_id = sparse_matrix_id,
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
                    {params.distance} {params.threads} {input.sparse_matrix_id} \
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
        sparse_matrix_id = sparse_matrix_id,
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
                    {input.sparse_matrix_id} \
                    > >(tee -a {log.stdout}) \
                    2> >(tee -a {log.stderr} >&2)
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
            exit 1
        fi
        """

rule filter_group_file:
    input:
        group = lambda wildcards: group_files
    output:
        "run_files/{gene}_group_file.txt"
    log:
        stdout="logs/filter_group_file/{gene}.out",
        stderr="logs/filter_group_file/{gene}.err"
    shell:
        """
        bash scripts/filter_group_file.sh {wildcards.gene} {output} \
            {input.group} \
            > >(tee -a {log.stdout}) \
            2> >(tee -a {log.stderr} >&2)
        """

# Format-specific stepwise conditional rules
rule spa_tests_stepwise_conditional_vcf:
    input:
        vcf="run_files/{gene}_{distance}_{maf}.vcf.bgz",
        vcf_csi="run_files/{gene}_{distance}_{maf}.vcf.bgz.csi",
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
        use_null_var_ratio=config["use_null_var_ratio"],
        P_T=config["conditioning_pvalue"]
    log:
        stdout="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.out",
        stderr="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.err"
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for vcf in {input.vcf}; do
            conda run --no-capture-output -n RSAIGE_vcf_version \
                bash scripts/stepwise_conditional_SAIGE.sh \
                $vcf {output} {input.model_file} {input.variance_file} \
                {input.sparse_matrix} $chr {params.use_null_var_ratio} \
                {params.P_T} \
                > >(tee -a {log.stdout}) \
                2> >(tee -a {log.stderr} >&2)
        done
        """

rule spa_tests_stepwise_conditional_plink:
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
        sparse_matrix_id=sparse_matrix_id,
        group_file="run_files/{gene}_group_file.txt",
    output:
        "run_files/{gene}_{trait}_{distance}_{maf}_string.txt"
    params:
        maf_common="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"],
        P_T=config["conditioning_pvalue"]
    log:
        stdout="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.out",
        stderr="logs/spa_tests_stepwise_conditional/{gene}_{trait}_{distance}_{maf}.err"
    shell:
        """
        set -euo pipefail
        chr=$(python scripts/extract_chromosome.py --ensembl_id \"{wildcards.gene}\")
        for plink_bed in {input.plink_bed}; do
            plink_fileset=$(echo "$plink_bed" | sed 's/\\.bed$//')
            conda run --no-capture-output -n RSAIGE_vcf_version \
                bash scripts/stepwise_conditional_SAIGE_plink.sh \
                $plink_fileset {output} {input.model_file} {input.variance_file} \
                {input.sparse_matrix} $chr {params.use_null_var_ratio} \
                {params.P_T} \
                > >(tee -a {log.stdout}) \
                2> >(tee -a {log.stderr} >&2)
        done
        """

# Format-specific conditional analysis rules
rule spa_tests_conditional_vcf:
    input:
        vcf = "run_files/{gene}_{distance}.vcf.bgz",
        vcf_csi = "run_files/{gene}_{distance}.vcf.bgz.csi",
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
    log:
        stdout="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.out",
        stderr="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.err"
    threads: 1
    shell:
        """
        set -euo pipefail
        conda run --no-capture-output -n RSAIGE_vcf_version \
            bash scripts/saige_step2_conditioning_check.sh \
            {input.vcf} {output} {params.min_mac} {input.model_file} \
            {input.variance_file} {input.sparse_matrix} {input.group_file} \
            {params.annotations_to_include} {input.conditioning_variants} \
            {params.max_MAF} {params.use_null_var_ratio} \
            > >(tee -a {log.stdout}) \
            2> >(tee -a {log.stderr} >&2)
        """

rule spa_tests_conditional_plink:
    input:
        plink_bim = "run_files/{gene}_{distance}.bim",
        plink_bed = "run_files/{gene}_{distance}.bed",
        plink_fam = "run_files/{gene}_{distance}.fam",
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
    log:
        stdout="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.out",
        stderr="logs/spa_tests_conditional/{gene}_{trait}_{distance}_{maf}.err"
    threads: 1
    shell:
        """
        set -euo pipefail
        plink_fileset=$(echo {input.plink_bed} | sed 's/\\.bed$//')
        conda run --no-capture-output -n RSAIGE_vcf_version \
            bash scripts/saige_step2_conditioning_check_plink.sh \
            $plink_fileset {output} {params.min_mac} {input.model_file} \
            {input.variance_file} {input.sparse_matrix} {input.group_file} \
            {params.annotations_to_include} {input.conditioning_variants} \
            {params.max_MAF} {params.use_null_var_ratio} \
            > >(tee -a {log.stdout}) \
            2> >(tee -a {log.stderr} >&2)
        """

rule combine_results:
    input:
        expand("saige_outputs/{gene_trait}_{distance}_saige_results_{maf}.txt",
               gene_trait=valid_gene_trait_pairs,
               distance=config["distance"],
               maf=config["maf"]),
    output:
        "brava_stepwise_conditional_analysis_results.txt"
    log:
        stdout="logs/combine_results/final_output.out",
        stderr="logs/combine_results/final_output.err"
    shell:
        """
        set -euo pipefail
        python scripts/combine_saige_outputs.py --out {output} \
        > >(tee -a {log.stdout}) \
        2> >(tee -a {log.stderr} >&2)
        """

# Rule order to ensure proper execution based on input format
if input_format == "vcf":
    ruleorder: spa_tests_stepwise_conditional_vcf > spa_tests_stepwise_conditional_plink
    ruleorder: spa_tests_conditional_vcf > spa_tests_conditional_plink
    ruleorder: filter_to_coding_gene_vcf > filter_to_coding_gene_plink
    ruleorder: filter_to_gene_vcf > filter_to_gene_plink
else:  # plink
    ruleorder: spa_tests_stepwise_conditional_plink > spa_tests_stepwise_conditional_vcf
    ruleorder: spa_tests_conditional_plink > spa_tests_conditional_vcf
    ruleorder: filter_to_coding_gene_plink > filter_to_coding_gene_vcf
    ruleorder: filter_to_gene_plink > filter_to_gene_vcf

ruleorder: identify_gene_start_stop > filter_to_coding_gene_plink
ruleorder: identify_gene_start_stop > filter_to_coding_gene_vcf
ruleorder: spa_tests_stepwise_conditional_vcf > filter_to_gene_vcf
ruleorder: spa_tests_stepwise_conditional_plink > filter_to_gene_plink
