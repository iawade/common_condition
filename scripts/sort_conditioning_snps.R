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
sort_conditoioning_snps <- function(snps) {
  cat("Sorting these condition SNPs:", snps, "\n")
  # ... do something more useful here ...
}

# Call the function with the argument from the command line
if (is.null(opt$condition)) {
  stop("Please provide a --condition argument.", call. = FALSE)
} else {
  sort_conditoioning_snps(opt$condition)
}
