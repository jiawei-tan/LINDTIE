################################################################################
# Module      : LINDTIE_compare_transcript
# Description : Performs differential expression analysis
# Copyright   : (c) Jia Wei Tan, Dec 2025
# License     : MIT
# Maintainer  : https://github.com/jiawei-tan
# Portability : POSIX
#
# Usage:
# Rscript LINDTIE_compare_transcript.R <transcript_counts_matrix_file> \
#   <oarfish_output_dir> \
#   <tx_ref_fasta> \
#   <outfile> \
#   --FDR=<value> --minCPM=<value> --minLogFC=<value>
#
# Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
################################################################################

library(dplyr)
library(data.table)

options(error = function() {
  logcat("An error occurred. Traceback has been captured below:\n")
  tb <- tryCatch(
    capture.output(traceback(2)),
    error = function(e) "No traceback available."
  )
  logcat(paste(tb, collapse = "\n"), "\n")
  q("no", status = 1, runLast = FALSE)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
if (length(file_arg) == 0) {
  stop(paste(
    "Unable to locate the script file path.",
    "Ensure the script is run with Rscript."
  ))
}
incl_path <- gsub(
  "--file=(.*)LINDTIE_compare_transcript.R",
  "\\1LINDTIE_de_methods.R",
  file_arg
)
source(incl_path, chdir = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Invalid arguments.\n\n")
  args <- c("--help")
}

if ("--help" %in% args) {
  cat("\nUsage:\n  Rscript LINDTIE_compare_transcript.R",
      "<case_name> <oarfish_output_dir> <tx_ref_fasta> <output>",
      "--FDR=<value> --minCPM=<value> --minLogFC=<value>\n\n")
  q(save = "no")
}

outfile <- args[4]
logfile <- file.path(dirname(outfile), "DE.log")

logcat <- function(...) {
  cat(..., file = logfile, append = TRUE)
  cat(...)
}

logcat("Starting LINDTIE_compare_transcript.R script...\n")
logcat("Arguments received: ", paste(args, collapse = " "), "\n")

# Set default values
FDR <- 0.05
minCPM <- 0.1
minLogFC <- 2

set_arg <- function(argname) {
  flag <- paste0("--", argname, "=")
  arg <- grep(flag, args, value = TRUE)
  if (length(arg) != 0) {
    value <- as.numeric(strsplit(arg, "=")[[1]][2])
    if (is.na(value)) {
      stop(paste("Invalid", argname, "value."))
    }
    assign(argname, value, envir = .GlobalEnv)
    logcat("Set ", argname, " to ", value, "\n")
  }
}

set_arg("FDR")
set_arg("minCPM")
set_arg("minLogFC")

case_name <- args[1]
oarfish_output_dir <- args[2]
tx_ref_fasta <- args[3]

logcat("Case name: ", case_name, "\n")
logcat("Oarfish output directory: ", oarfish_output_dir, "\n")
logcat("Transcript FASTA: ", tx_ref_fasta, "\n")
logcat("Output file: ", outfile, "\n")

logcat("Running differential expression analysis...\n")
de_results <- run_edgeR(case_name, oarfish_output_dir,
                        outdir = dirname(outfile), cpm_cutoff = minCPM,
                        qval = FDR, min_logfc = minLogFC)

write.table(de_results, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
logcat("Final results written to ", outfile, "\n")
logcat("LINDTIE_compare_transcript.R script completed successfully.\n")
