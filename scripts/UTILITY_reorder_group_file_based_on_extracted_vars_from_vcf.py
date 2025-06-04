import gzip

# --- Setup paths
input_file = "../groups/combined_group_file.txt.gz"
output_file = "filtered_group_file.txt.gz"
id_file = "all_var_ids.txt"

# --- Load valid variant IDs
print(f"[INFO] Loading valid variant IDs from: {id_file}")
with open(id_file) as f:
    valid_ids = set(line.strip() for line in f)
print(f"[INFO] Loaded {len(valid_ids):,} valid variant IDs.")

# --- Open input/output files
print(f"[INFO] Processing group file: {input_file}")
print(f"[INFO] Writing filtered output to: {output_file}")
genes_kept = 0
genes_skipped = 0

with gzip.open(input_file, "rt") as fin, gzip.open(output_file, "wt") as fout:
    while True:
        var_line = fin.readline()
        anno_line = fin.readline()

        if not var_line or not anno_line:
            break  # End of file

        var_fields = var_line.strip().split()
        anno_fields = anno_line.strip().split()

        gene_id = var_fields[0]
        variants = var_fields[2:]
        annotations = anno_fields[2:]

        # Filter by valid variant IDs
        filtered = [(v, a) for v, a in zip(variants, annotations) if v in valid_ids]

        if filtered:
            kept_vars, kept_annos = zip(*filtered)
            fout.write(f"{gene_id} var {' '.join(kept_vars)}\n")
            fout.write(f"{gene_id} anno {' '.join(kept_annos)}\n")
            genes_kept += 1
        else:
            genes_skipped += 1

print(f"[DONE] Finished filtering.")
print(f"[STATS] Genes written: {genes_kept:,}")
print(f"[STATS] Genes skipped (no valid variants): {genes_skipped:,}")
