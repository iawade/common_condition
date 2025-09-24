# Snakemake Pipeline for BRaVa Pilot Study Common Variant Conditional Analysis
# Supports both VCF and PLINK input formats

configfile: "config.yaml"

input_format = config["input_format"]  # "vcf" or "plink"
sparse_matrix = config["sparse_matrix"]
sparse_matrix_id = f"{sparse_matrix}.sampleIDs.txt"
gene_trait_pairs_to_test_with_conditioning_variants = config["gene_trait_pairs_to_test_with_conditioning_variants"]
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
print(f"Input files Sample: {input_files[:5]}")
print(f"Group files Sample: {group_files[:5]}")

distance = config["distance"]
min_mac = config["min_mac"]
annotations_to_include = config["annotations_to_include"]

import pandas as pd
import json
import re
from scripts.extract_chromosome import get_gene_chr

# Extract chromosomes based on input format
pattern = r'\.(chr[0-9X]+)\.'
if input_format == "plink":
    matches = [re.search(pattern, s) for s in plink_bim_files]
else:  # vcf
    matches = [re.search(pattern, s) for s in input_files]
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
    trait_in_variance = any(re.search(pattern, vf) for vf in variance_files)
    if trait_in_model and trait_in_variance:
        available_traits.add(pid)

conditioning_jobs = [job for job in conditioning_jobs if job['Trait'] in available_traits]

for job in conditioning_jobs:
    filename = f"final_run_files/{job['Gene']}_{job['Trait']}_{job['MAF_cutoff_for_conditioning_variants']}_extract.txt"
    with open(filename, "w") as f:
        f.writelines(f"{v}\n" for v in job['conditioning_variants'].split(","))
    filename = f"final_run_files/{job['Gene']}_{job['Trait']}_{job['MAF_cutoff_for_conditioning_variants']}_extract.bed"
    with open(filename, "w") as f:
        variants = job['conditioning_variants'].split(",")
        for v in variants:
            chrom, pos, ref, alt = v.split(":")
            _ = f.write(f"{chrom}\t{pos}\t{pos}\n")

genes = sorted(set(job['Gene'] for job in conditioning_jobs))

# Debugging: Print valid gene-trait pairs
print(f"Filtered {len(conditioning_jobs)} conditioning jobs with available model/variance files.")

# Format-specific output files for rule all
def get_format_outputs():
    base_outputs = [
        expand([
            "final_run_files/{gene}_{trait}_{maf}_extract.txt",
            "final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string.txt",
            "final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}_{distance}.txt"],
        zip,
        gene=[job['Gene'] for job in conditioning_jobs],
        trait=[job['Trait'] for job in conditioning_jobs],
        maf=[job['MAF_cutoff_for_conditioning_variants'] for job in conditioning_jobs],
        distance=[config["distance"]] * len(conditioning_jobs)
        ),
        expand("final_run_files/{gene}_group_file.txt", gene=genes),
        expand("final_run_files/bed/{gene}.bed", gene=genes),
        expand("final_run_files/bed/expanded_regions_{gene}.bed", gene=genes),
        "brava_final_conditional_analysis_results.txt"
    ]
    
    # Add format-specific outputs
    if input_format == "plink":
        base_outputs.extend([
            expand("final_run_files/{gene}_{distance}.bim", 
                gene=genes, distance=config["distance"]),
            expand("final_run_files/{gene}_{distance}.bed", 
                gene=genes, distance=config["distance"]),
            expand("final_run_files/{gene}_{distance}.fam", 
                gene=genes, distance=config["distance"])
        ])
    else:  # vcf
        base_outputs.extend([
            expand("final_run_files/{gene}_{distance}.vcf.bgz", 
                gene=genes, distance=config["distance"]),
            expand("final_run_files/{gene}_{distance}.vcf.bgz.csi", 
                gene=genes, distance=config["distance"])
        ])
    
    return base_outputs

# Target Rule for Completion of Pipeline
rule all:
    input:
        get_format_outputs()

rule filter_group_file:
    input:
        group = lambda wildcards: group_files
    output:
        "final_run_files/{gene}_group_file.txt"
    log:
        stdout="logs/final_filter_group_file/{gene}.out",
        stderr="logs/final_filter_group_file/{gene}.err"
    shell:
        """
        bash scripts/filter_group_file.sh {wildcards.gene} {output} {input.group} > {log.stdout} 2> {log.stderr}
        """

rule identify_gene_start_stop:
    output:
        "final_run_files/bed/{gene}.bed",
        "final_run_files/bed/expanded_regions_{gene}.bed"
    params:
        distance=distance,
        outfolder="final_run_files/bed"
    log:
        stdout="logs/identify_gene_start_stop/{gene}.out",
        stderr="logs/identify_gene_start_stop/{gene}.err"
    shell:
        """
        set -euo pipefail
        python scripts/start_end_query.py --ensembl_id \"{wildcards.gene}\" --folder {params.outfolder} > {log.stdout} 2> {log.stderr}
        bash scripts/expand_coding_region.sh {wildcards.gene} {params.distance} {params.outfolder} >> {log.stdout} 2>> {log.stderr}
        """

# VCF-specific rules
rule filter_to_gene_vcf:
    input:
        vcf = lambda wildcards: input_files if input_format == "vcf" else [],
        bed = "final_run_files/bed/expanded_regions_{gene}.bed"
    output:
        "final_run_files/{gene}_{distance}.vcf.bgz",
        "final_run_files/{gene}_{distance}.vcf.bgz.csi"
    params:
        distance=distance,
        threads=config["threads"],
        outfolder="final_run_files"
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
                bash scripts/filter_to_gene_vcf.sh $vcf {wildcards.gene} {params.distance} {params.threads} >> {log.stdout} 2>> {log.stderr}
            fi
        done

        if [[ -z "$matched_vcf" ]]; then
            echo "ERROR: No matching VCF found for chromosome $chr"
        fi
        """

# PLINK-specific rules
rule filter_to_gene_plink:
    input:
        plink_bim = lambda wildcards: plink_bim_files if input_format == "plink" else [],
        plink_bed = lambda wildcards: plink_bed_files if input_format == "plink" else [],
        plink_fam = lambda wildcards: plink_fam_files if input_format == "plink" else [],
        sparse_matrix_id = sparse_matrix_id,
        regions = "final_run_files/bed/expanded_regions_{gene}.bed"
    output:
        bim = "final_run_files/{gene}_{distance}.bim",
        bed = "final_run_files/{gene}_{distance}.bed",
        fam = "final_run_files/{gene}_{distance}.fam"
    params:
        distance=distance,
        threads=config["threads"],
        outfolder="final_run_files"
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
                bash scripts/filter_to_gene_plink.sh $plink_fileset {wildcards.gene} {params.distance} {params.threads} {input.sparse_matrix_id} >> {log.stdout} 2>> {log.stderr}
            fi
        done

        if [[ -z "$matched_plink" ]]; then
            echo "ERROR: No matching plink fileset (.bim/.bed/.fam) found for chromosome $chr"
        fi
        """

# Format-specific pruning rules
rule prune_to_independent_conditioning_variants_vcf:
    input:
        vcf = "final_run_files/{gene}_{distance}.vcf.bgz",
        vcf_csi = "final_run_files/{gene}_{distance}.vcf.bgz.csi",
        conditioning_variants = "final_run_files/{gene}_{trait}_{maf}_extract.txt",
        conditioning_variants_bed = "final_run_files/{gene}_{trait}_{maf}_extract.bed",
        group_file="final_run_files/{gene}_group_file.txt"
    output:
        "final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string.txt" 
    params:
        file="final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string"
    shell:
        """
        set -euo pipefail

        TMPFILE=$(mktemp)
        # Define a bed file first
        plink2 --vcf {input.vcf} --extract range {input.conditioning_variants_bed} \
          --make-bed --out ${{TMPFILE}} || true
        if [[ -f ${{TMPFILE}}.bim ]]; then

            # Sort the .bim file variant ID
            awk 'BEGIN{{OFS="\t"}} {{ $2 = "chr" $1 ":" $4 ":" $6 ":" $5; print }}' ${{TMPFILE}}.bim > ${{TMPFILE}}.bim.tmp
            mv ${{TMPFILE}}.bim.tmp ${{TMPFILE}}.bim

            # Extract the exact variants
            plink2 --bfile ${{TMPFILE}} \
                --extract {input.conditioning_variants} \
                --make-bed --out ${{TMPFILE}}.tmp || true

            nvar=$(wc -l < ${{TMPFILE}}.tmp.bim)

            if [[ $nvar -eq 1 ]]; then
                cut -f2 ${{TMPFILE}}.tmp.bim > {output}
            else
                plink2 --bfile ${{TMPFILE}}.tmp \
                  --indep-pairwise 50 5 0.9 \
                  --out {params.file} || true
                # Finally, create a comma separated string from this
                if [[ -f {params.file}.prune.in ]]; then
                    paste -sd, {params.file}.prune.in > {output}
                fi
            fi
        else
            # No conditioning variants present in the vcf file
            touch {output}
        fi
        """

rule prune_to_independent_conditioning_variants_plink:
    input:
        plink_bim = "final_run_files/{gene}_{distance}.bim",
        plink_bed = "final_run_files/{gene}_{distance}.bed",
        plink_fam = "final_run_files/{gene}_{distance}.fam",
        conditioning_variants = "final_run_files/{gene}_{trait}_{maf}_extract.txt",
        group_file="final_run_files/{gene}_group_file.txt"
    output:
        "final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string.txt"
    params:
        file="final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string"
    shell:
        """
        set -euo pipefail

        TMPFILE=$(mktemp)
        plink_fileset=$(echo "{input.plink_bed}" | sed 's/\\.bed$//')
        plink2 --bfile $plink_fileset \
            --extract {input.conditioning_variants} \
            --make-bed --out ${{TMPFILE}} || true
        
        if [[ -f ${{TMPFILE}}.bim ]]; then
            nvar=$(wc -l < ${{TMPFILE}}.bim)
            if [[ $nvar -eq 1 ]]; then
                cut -f2 ${{TMPFILE}}.bim > {output}
            else
                plink2 --bfile ${{TMPFILE}} \
                  --indep-pairwise 50 5 0.9 \
                  --out {params.file} || true
                # Finally, create a comma separated string from this
                if [[ -f {params.file}.prune.in ]]; then
                    paste -sd, {params.file}.prune.in > {output}
                fi
            fi
        else
            # No variants to prune, create an empty output
            touch {output}
        fi
        """

# Format-specific conditional analysis rules
rule spa_tests_conditional_vcf:
    input:
        vcf = "final_run_files/{gene}_{distance}.vcf.bgz",
        vcf_csi = "final_run_files/{gene}_{distance}.vcf.bgz.csi",
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
        conditioning_variants="final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string.txt"
    output:
        "final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}_{distance}.txt"
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    log:
        stdout="logs/final_spa_tests_conditional/{gene}_{trait}_{maf}_{distance}.out",
        stderr="logs/final_spa_tests_conditional/{gene}_{trait}_{maf}_{distance}.err"
    threads: 4
    shell:
        """
        set -euo pipefail
        conda run --no-capture-output -n RSAIGE_vcf_version \
            bash scripts/saige_step2_conditioning_check.sh \
            {input.vcf} {output} {params.min_mac} {input.model_file} \
            {input.variance_file} {input.sparse_matrix} {input.group_file} \
            {params.annotations_to_include} {input.conditioning_variants} \
            {params.max_MAF} {params.use_null_var_ratio} > {log.stdout} 2> {log.stderr}
        """

rule spa_tests_conditional_plink:
    input:
        plink_bim = "final_run_files/{gene}_{distance}.bim",
        plink_bed = "final_run_files/{gene}_{distance}.bed",
        plink_fam = "final_run_files/{gene}_{distance}.fam",
        model_file = lambda wildcards: [
            mf for mf in model_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', mf)
        ],
        variance_file=lambda wildcards: [
            vf for vf in variance_files
            if re.search(rf'(?:^|[/_.\-]){re.escape(wildcards.trait)}(?=[/_.\-])', vf)
        ],
        sparse_matrix=sparse_matrix,
        group_file="final_run_files/{gene}_group_file.txt",
        conditioning_variants="final_run_files/{gene}_{trait}_{maf}_{distance}_ld_pruned_string.txt"
    output:
        "final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}_{distance}.txt" 
    params:
        min_mac=min_mac,
        annotations_to_include=annotations_to_include,
        max_MAF="{maf}",
        use_null_var_ratio=config["use_null_var_ratio"]
    log:
        stdout="logs/final_spa_tests_conditional/{gene}_{trait}_{maf}_{distance}.out",
        stderr="logs/final_spa_tests_conditional/{gene}_{trait}_{maf}_{distance}.err"
    threads: 4
    shell:
        """
        set -euo pipefail
        plink_fileset=$(echo {input.plink_bed} | sed 's/\\.bed$//')
        conda run --no-capture-output -n RSAIGE_vcf_version \
            bash scripts/saige_step2_conditioning_check_plink.sh \
            $plink_fileset {output} {params.min_mac} {input.model_file} \
            {input.variance_file} {input.sparse_matrix} {input.group_file} \
            {params.annotations_to_include} {input.conditioning_variants} \
            {params.max_MAF} {params.use_null_var_ratio} > {log.stdout} 2> {log.stderr}
        """

rule combine_results:
    input:
        expand("final_saige_outputs/{gene}_{trait}_saige_conditioned_results_{maf}_{distance}.txt",
        zip,
        gene=[job['Gene'] for job in conditioning_jobs],
        trait=[job['Trait'] for job in conditioning_jobs],
        maf=[job['MAF_cutoff_for_conditioning_variants'] for job in conditioning_jobs],
        distance=[config["distance"]] * len(conditioning_jobs)
        ),
    output:
        "brava_final_conditional_analysis_results.txt"
    log:
        stdout="logs/final_combine_results/final_output.out",
        stderr="logs/final_combine_results/final_output.err"
    shell:
        """
        set -euo pipefail
        python scripts/combine_saige_outputs.py --final --out {output} > {log.stdout} 2> {log.stderr}
        """

# Rule order to ensure proper execution based on input format
if input_format == "vcf":
    ruleorder: prune_to_independent_conditioning_variants_vcf > prune_to_independent_conditioning_variants_plink
    ruleorder: spa_tests_conditional_vcf > spa_tests_conditional_plink
    ruleorder: filter_to_gene_vcf > filter_to_gene_plink
else:  # plink
    ruleorder: prune_to_independent_conditioning_variants_plink > prune_to_independent_conditioning_variants_vcf
    ruleorder: spa_tests_conditional_plink > spa_tests_conditional_vcf
    ruleorder: filter_to_gene_plink > filter_to_gene_vcf

ruleorder: identify_gene_start_stop > filter_to_gene_plink
ruleorder: identify_gene_start_stop > filter_to_gene_vcf