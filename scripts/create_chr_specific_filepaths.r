#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages(library(optparse))

# Define command line options
option_list <- list(
  make_option(c("-c", "--chr"), type = "character", default = NULL,
              help = "chromosome to restrict to", metavar = "character"),
  make_option(c("-p", "--pvalue_T"), type = "numeric", default = NULL,
              help = "alter the p-value threshold in the config")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(data.table)
library(dplyr)
library(yaml)

config <- read_yaml("config.yaml")

if (is.null(opt$chr)) {
  cat("No chromosome restriction")
} else {
  opt$chr <- ifelse(grepl("^chr", opt$chr), opt$chr, paste0("chr", opt$chr))
  if (!(opt$chr %in% paste0("chr", c(seq(1,22), "X")))) {
    stop("Please provide a --chr argument in 1..22, X.", call. = FALSE)
  }

  # Determine the vcf/plink filepaths
  paths_to_change <- c("list_of_vcf_files", "list_of_plink_files",
    "list_of_group_files")

  for (path in paths_to_change) {
    if (path %in% names(config)) {
      if (!is.null(config[[path]])) {
        files <- fread(config[[path]], header=FALSE)$V1
        match <- grep(paste0("(^|[._-])", opt$chr, "($|[._-])"), files, perl = TRUE)
        if (length(match) != 1) {
          stop(paste("Multiple or no chr matches in the", path, "paths"),
            call. = FALSE)
        }
        out <- paste0(gsub("list_of_", "", path), ".txt")
        fwrite(data.table(files=files[match]), quote=FALSE, col.names=FALSE,
          file=out)
        config[[path]] <- out
      }
    }
  }
}

if (is.null(opt$p)) {
  cat("No alteration to the passed p-value threshold")
} else {
  config[["conditioning_pvalue"]] <- opt$p
}

# Rewrite the config
write_yaml(config, "config.yaml")
