#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages(library(optparse))

# Define command line options
option_list <- list(
  make_option(c("-c", "--condition"), type = "character", default = NULL,
              help = "Condition name to use", metavar = "character")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Define a function that does something with the condition
sort_conditioning_snps <- function(snps)
{
  message("Sorting these conditioning SNPs: ", snps, "\n")
  snps_vec <- strsplit(snps, split=",")[[1]]
  if (length(snps_vec) > 1) {
    return(paste(snps_vec[order(as.integer(unlist(lapply(strsplit(snps_vec, split=":"), `[`, 2))))], collapse=","))
  } else {
    return(snps)
  }
}

# Call the function with the argument from the command line
if (is.null(opt$condition)) {
  stop("Please provide a --condition argument.", call. = FALSE)
} else {
  snps <- sort_conditioning_snps(opt$condition)
  cat(snps)
}
