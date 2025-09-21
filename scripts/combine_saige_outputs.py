import pandas as pd
import glob
import re
import os
import argparse

def combine_outputs(out, final=False):

    # Define the file pattern
    if (final):
        file_pattern = "final_saige_outputs/*_*_saige_*results_*.txt"
        file_single_pattern = "final_saige_outputs/*_*_saige_*results_*.singleAssoc.txt"
        file_string_pattern = "final_run_files/*_*_*_ld_pruned_string.txt"
    else:
        file_pattern = "saige_outputs/*_*_*_saige_results_*.txt"
        file_single_pattern = "saige_outputs/*_*_*_saige_results_*.singleAssoc.txt"
        file_string_pattern = "run_files/*_*_*_*_string.txt"
    
    # SAIGE outputs
    files = [f for f in glob.glob(file_pattern)]
    files_single = [f for f in glob.glob(file_single_pattern)]
    files = [item for item in files if item not in files_single]

    # Conditioning variants
    files_string = [f for f in glob.glob(file_string_pattern)]

    dfs = []
    df_singles = []
    df_conditionings = []

    # Process each file
    for file in files:
        # Extract ancestry, trait, variant class and mode from filename
        if (final):
            match = re.search(r"final_saige_outputs/(.*?)_(.*?)_saige_conditioned_results_(.*?)\.txt", file)
        else:
            match = re.search(r"saige_outputs/(.*?)_(.*?)_.*?_saige_results_(.*?)\.txt", file)

        if not match:
            continue
        if os.stat(file).st_size == 0:
            continue
        gene, trait, maf_cutoff = match.groups()

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
        print(f"results file is not empty: saved as {out}")
    else:
        open(out, 'a').close()
        print(f"results file is empty: saved as {out}")

    # Process each file
    for file in files_single:
        # Extract ancestry, trait, variant class and mode from filename
        if (final):
            match = re.search(r"final_saige_outputs/(.*?)_(.*?)_saige_conditioned_results_(.*?)\.txt\.singleAssoc\.txt", file)
        else:
            match = re.search(r"saige_outputs/(.*?)_(.*?)_.*?_saige_results_(.*?)\.txt\.singleAssoc\.txt", file)

        if not match:
            continue
        if os.stat(file).st_size == 0:
            continue
        gene, trait, maf_cutoff = match.groups()

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
        combined_df.to_csv(out + ".singleAssoc.txt", sep="\t", index=False)
        print(f"singleAssoc file is not empty: saved as {out}.singleAssoc.txt")
    else:
        open(out + ".singleAssoc.txt", 'a').close()
        print(f"singleAssoc file is empty: saved as {out}.singleAssoc.txt")

    # Process each file
    for file in files_string:
        # Extract ancestry, trait, variant class and mode from filename
        if (final):
            match = re.search(r"final_run_files/(.*?)_(.*?)_(.*?)_ld_pruned_string\.txt", file)
        else:
            match = re.search(r"run_files/(.*?)_(.*?)_.*?_(.*?)_string\.txt", file)

        if not match:
            continue
        if os.stat(file).st_size == 0:
            continue
        gene, trait, maf_cutoff = match.groups()

        if os.path.getsize(file) == 0:
            print("File is empty")
        else:
            # Read file
            df = pd.read_csv(file, sep="\t", header=None, names=["conditioning_variants"])
            # Add extracted metadata as new columns
            df["Gene"] = gene
            df["Trait"] = trait
            df["MAF_cutoff_for_conditioning_variants"] = maf_cutoff
            df_conditionings.append(df)

    if (len(df_conditionings) > 0):
        # Combine all dataframes
        combined_df = pd.concat(df_conditionings, ignore_index=True)
        # Save to a new file
        combined_df.to_csv(out + ".conditioning.variants.txt", sep="\t", index=False)
        print(f"conditioning variants file is not empty: saved as {out}.conditioning.variants.txt")
    else:
        open(out + ".conditioning.variants.txt", 'a').close()
        print(f"conditioning variants file is empty: saved as {out}.conditioning.variants.txt")

    print(f"Files merged and sorted. Output saved as {out}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Combine all of the output files together.')
    parser.add_argument('--out', required=True, help='Name of the output file.')
    args = parser.parse_args()

    # Run the function with provided arguments
    combine_outputs(args.out)