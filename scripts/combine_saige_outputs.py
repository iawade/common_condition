import pandas as pd
import glob
import re
import os
import argparse

def combine_outputs(out):

    # Define the file pattern
    file_pattern = "saige_outputs/*_*_saige_results_*.txt"
    files = [f for f in glob.glob(file_pattern)]
    file_single_pattern = "saige_outputs/*_*_saige_results_*.singleAssoc.txt"
    files_single = [f for f in glob.glob(file_single_pattern)]
    files = [item for item in files if item not in files_single]

    dfs = []
    df_singles = []

    # Process each file
    for file in files:
        # Extract ancestry, trait, variant class and mode from filename
        match = re.search(r"saige_outputs/(.*?)_(.*?)_(.*?)_saige_results_(.*?)\.txt", file)

        if not match:
            continue
        if os.stat(file).st_size == 0:
            continue
        gene, trait, distance, maf_cutoff = match.groups()

        if os.path.getsize(file) == 0:
            print("File is empty")
        else:
            # Read file
            df = pd.read_csv(file, sep="\t")
            # Add extracted metadata as new columns
            df["Gene_sanity_check"] = gene
            df["Trait"] = trait
            df["MAF_cutoff_for_conditioning_variants"] = maf_cutoff
            dfs.append(df)

    if (len(dfs) > 0):
        # Combine all dataframes
        combined_df = pd.concat(dfs, ignore_index=True)
        # Save to a new file
        combined_df.to_csv(out, sep="\t", index=False)
    else:
        open(out, 'a').close()

    # Process each file
    for file in files_single:
        # Extract ancestry, trait, variant class and mode from filename
        match = re.search(r"saige_outputs/(.*?)_(.*?)_(.*?)_saige_results_(.*?)\.singleAssoc\.txt", file)

        if not match:
            continue
        if os.stat(file).st_size == 0:
            continue
        gene, trait, distance, maf_cutoff = match.groups()

        if os.path.getsize(file) == 0:
            print("File is empty")
        else:
            # Read file
            df = pd.read_csv(file, sep="\t")
            # Add extracted metadata as new columns
            df["Gene_sanity_check"] = gene
            df["Trait"] = trait
            df["MAF_cutoff_for_conditioning_variants"] = maf_cutoff
            df_singles.append(df)

    if (len(df_singles) > 0):
        # Combine all dataframes
        combined_df = pd.concat(df_singles, ignore_index=True)
        # Save to a new file
        combined_df.to_csv(out + "singleAssoc.txt", sep="\t", index=False)
    else:
        open(out + "singleAssoc.txt", 'a').close()

    print(f"Files merged and sorted. Output saved as {out}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Combine all of the output files together.')
    parser.add_argument('--out', required=True, help='Name of the output file.')
    args = parser.parse_args()

    # Run the function with provided arguments
    combine_outputs(args.out)