#!/usr/bin/env Rscript

# Usage:
# Rscript filter_group_file.R <gene> <output_prefix> <annotation_file1[,annotation_file2,...]> <setlist_file1[,setlist_file2,...]>

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript filter_group_file.R <gene> <output_prefix> <annotation_file1[,annotation_file2,...]> <setlist_file1[,setlist_file2,...]>")
}

gene <- args[1]
outfile_prefix <- args[2]

# Split comma-separated lists
annotation_files <- strsplit(args[3], ",")[[1]]
setlist_files <- strsplit(args[4], ",")[[1]]

message("Filtering for gene: ", gene)
message("Annotation files: ", paste(annotation_files, collapse = ", "))
message("Setlist files: ", paste(setlist_files, collapse = ", "))

# Helper function to read a file safely
safe_fread <- function(file) {
  if (!file.exists(file)) {
    warning("File not found: ", file)
    return(NULL)
  }
  fread(file, sep = " ", header = FALSE, data.table = TRUE)
}

# Filter and combine multiple files
filter_files <- function(files, gene_col) {
  combined <- lapply(files, function(f) {
    dt <- safe_fread(f)
    if (!is.null(dt)) {
      filtered <- dt %>% filter(.data[[gene_col]] == gene)
      message(basename(f), ": matched ", nrow(filtered), " rows")
      return(filtered)
    }
    NULL
  })
  rbindlist(combined, fill = TRUE)
}

# V2 for annotation, V1 for setlist
dt_annot <- filter_files(annotation_files, "V2")
dt_setlist <- filter_files(setlist_files, "V1")

message("Total annotation rows matched: ", nrow(dt_annot))
message("Total setlist rows matched: ", nrow(dt_setlist))

# Write outputs
fwrite(dt_annot, file = paste0(outfile_prefix, ".annotation.txt"),
       sep = " ", quote = FALSE, col.names = FALSE)
fwrite(dt_setlist, file = paste0(outfile_prefix, ".setlist.txt"),
       sep = " ", quote = FALSE, col.names = FALSE)

message("✔ Done. Wrote filtered files with prefix: ", outfile_prefix)
