import pandas as pd
import glob
import re
import os

# Define the file pattern
file_pattern = "saige_outputs/*_*_saige_results_*.txt"
files = [f for f in glob.glob(file_pattern)]

dfs = []

# Process each file
for file in files:
    # Extract ancestry, trait, variant class and mode from filename
    match = re.search(r"saige_outputs/(.*?)_(.*?)_(.*?)_saige_results_(.*?)\.txt", file)
    if not match:
        continue
    if os.stat(file).st_size == 0:
        continue
    gene, trait, distance, maf_cutoff = match.groups()

    # Read file
    df = pd.read_csv(file, sep="\t")

    # Add extracted metadata as new columns
    df["Gene_sanity_check"] = gene
    df["Trait"] = trait
    df["MAF_cutoff_for_conditioning_variants"] = maf_cutoff

    dfs.append(df)

# Combine all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Save to a new file
combined_df.to_csv("brava_conditional_analysis_results.txt", sep="\t", index=False)

print("Files merged and sorted. Output saved as 'brava_conditional_analysis_results.txt'.")