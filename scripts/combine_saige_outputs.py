import pandas as pd
import glob
import re

# Define the file pattern
file_pattern = "saige_outputs/*_*_saige_results_*.txt"
files = [f for f in glob.glob(file_pattern)]

dfs = []

# Process each file
for file in files:
    # Extract ancestry, trait, variant class and mode from filename
    match = re.search(r"saige_outputs/(.*?)_(.*?)_saige_results_(.*?)\.txt", file)
    if not match:
        continue
    gene, trait, maf_cutoff = match.groups()

    # Read file
    df = pd.read_csv(file, sep="\t")

    # Add extracted metadata as new columns
    df["Gene"] = gene
    df["Trait"] = trait
    df["MAF_cutoff_for_common_variants"] = maf_cutoff

    # Apply filtering if columns exist
    df = df[df["AC_Allele2"] >= 5]

    if "N_case" in df.columns:
        df["N_case"] = pd.to_numeric(df["N_case"], errors="coerce")
        df = df[df["N_case"] >= 100]

    if "AF_Allele2" in df.columns:
        df["AF_Allele2"] = pd.to_numeric(df["AF_Allele2"], errors="coerce")
        df = df[df["AF_Allele2"] < 0.5]  # Keep only if AF_Allele2 < 0.5


    dfs.append(df)

# Combine all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Convert p.value to numeric (handle potential errors)
combined_df["p.value"] = pd.to_numeric(combined_df["p.value"], errors="coerce")

# Sort by p.value (ascending)
sorted_df = combined_df.sort_values(by="p.value", ascending=True)

# Save to a new file
sorted_df.to_csv("brava_conditional_analysis_results.txt", sep="\t", index=False)

print("Files merged and sorted. Output saved as 'brava_conditional_analysis_results.txt'.")