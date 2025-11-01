#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages(library(optparse))

# Define command line options
option_list <- list(
  make_option(c("-m", "--modelfile"), type = "character", default = NULL,
              help = "modelfile to use to filter samples to", metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "output file for column of sampleIDs", metavar = "character")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Define a function that does something with the condition
samples <- function(modelfile)
{
  message("Determining samples IDs from modelfile: ", modelfile, "\n")
  load(modelfile)
  return(modglmm$sampleID)
}

# Call the function with the argument from the command line
if (is.null(opt$condition)) {
  stop("Please provide a --modelfile argument.", call. = FALSE)
} else {
  samples <- samples(opt$modelfile)
  fwrite(samples, file=opt$outfile, sep="\t", quote=FALSE)
}
