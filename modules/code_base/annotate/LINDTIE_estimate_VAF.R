################################################################################
# Module      : LINDTIE_estimate_VAF
# Description : Estimates VAF and performs Fisher's Exact Test for Enrichment
# Copyright   : (c) Jia Wei Tan, Dec 2025
################################################################################

library(tximport)
library(data.table)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6) stop("Usage: Rscript LINDTIE_estimate_VAF.R <matrix> <quant> <cinfo> <fasta> <tx2gene> <outfile>")

transcript_counts_matrix_file <- args[1]
quant_file <- args[2]
contig_info_file <- args[3]
tx_ref_fasta <- args[4]
tx2gene_file <- args[5]
outfile <- args[6]

#############################################################
# Load data
#############################################################
print("Loading data...")
# Load abundance (TPM) and counts (Raw) from Oarfish
txi <- tximport(quant_file, type="oarfish", countsFromAbundance = "no", txOut=TRUE, dropInfReps=TRUE)
transcript_counts_matrix <- fread(transcript_counts_matrix_file)
cinfo <- fread(contig_info_file)
tx2g <- fread(tx2gene_file, col.names = c("transcript", "gene"))
txs <- fread(tx_ref_fasta, header = FALSE, sep = "\n")

#############################################################
# Prepare data
#############################################################
print("Preparing data...")
# Compute VAF for ALL annotated contigs (previously restricted to
# variant_of_interest == TRUE). VAF is gene-relative, so contigs without an
# overlapping gene are still excluded downstream by filter(gene != "").
# cinfo is used as-is; no variant_of_interest subsetting.

# Gene mapping logic
c2g <- data.frame(contig_id = cinfo$contig_id, gene=cinfo$overlapping_genes)
split_genes <- sapply(c2g$gene, function(x){strsplit(x, "\\||:")})
c2g <- data.frame(transcript = rep(c2g$contig_id, sapply(split_genes, length)),
                  gene = as.character(unlist(split_genes))) %>% distinct() %>% filter(gene != "")

colnames(transcript_counts_matrix)[colnames(transcript_counts_matrix) == "transcript_id"] <- "transcript"

# Reference cleaning logic
clean_tx_names <- function(names) {
    ifelse(names %like% "hg38_ncbiRefSeq", 
           sapply(gsub("hg38_ncbiRefSeq_", "", names), function(x) strsplit(x, "\\.")[[1]][1]), 
           names)
}

transcript_counts_matrix$transcript <- clean_tx_names(transcript_counts_matrix$transcript)
tx2g$transcript <- clean_tx_names(tx2g$transcript)
tx2g <- distinct(tx2g)

# Match non-novel contig transcripts to genes and combine with novel contig mappings
transcript2g <- right_join(transcript_counts_matrix, tx2g, by = "transcript")
tx2g <- transcript2g[, c("transcript", "gene")]
tx2g <- distinct(rbind(tx2g, c2g))

# Get reference list from fasta
txs_list <- clean_tx_names(gsub(">", "", sapply(txs$V1[grep("^>", txs$V1)], function(x) strsplit(x, " ")[[1]][1])))
tx_transcript <- transcript_counts_matrix
tx_transcript$novel <- !tx_transcript$transcript %in% txs_list
novel_contig_transcript <- tx_transcript[, all(novel), by = "transcript"]
novel_contig_transcript <- unique(novel_contig_transcript$transcript[novel_contig_transcript$V1])
novel_contigs <- unique(tx_transcript$transcript[tx_transcript$transcript %in% novel_contig_transcript])

# Identify Case vs Control columns
matrix_cols <- setdiff(colnames(transcript_counts_matrix), "transcript")
case_sample_id <- sub("\\.quant(\\.gz)?$", "", basename(quant_file))
if (case_sample_id %in% matrix_cols) {
    control_sample_cols <- setdiff(matrix_cols, case_sample_id)
} else {
    control_sample_cols <- matrix_cols
    if (length(control_sample_cols) > 0) {
        message("WARNING: Case sample id not found in transcript count matrix columns; treating all samples as controls.")
    }
}

# Calculate Control counts
if (length(control_sample_cols) > 0) {
    control_counts_summary <- transcript_counts_matrix %>%
        select(transcript, all_of(control_sample_cols)) %>%
        mutate(control_contig_count = rowSums(across(all_of(control_sample_cols)), na.rm = TRUE)) %>%
        select(transcript, control_contig_count)
} else {
    control_counts_summary <- transcript_counts_matrix %>%
        select(transcript) %>%
        mutate(control_contig_count = 0)
}

control_gene_counts <- left_join(control_counts_summary, tx2g, by = "transcript")
control_gene_total <- control_gene_counts %>%
    group_by(gene) %>%
    summarise(control_gene_count_total = sum(control_contig_count, na.rm = TRUE), .groups = "drop")
control_gene_wt <- control_gene_counts %>%
    filter(!transcript %in% novel_contigs) %>%
    group_by(gene) %>%
    summarise(control_gene_count_wt = sum(control_contig_count, na.rm = TRUE), .groups = "drop")

#############################################################
# Extract Raw Case Counts
#############################################################
print("Extracting raw case counts...")
case_counts_raw <- data.frame(
    num_reads_case = as.numeric(txi$counts[, 1]), 
    transcript = clean_tx_names(rownames(txi$counts))
)

case_wt_raw <- inner_join(case_counts_raw, tx2g, by = "transcript")
case_wt_gene_sum <- case_wt_raw %>%
    filter(!transcript %in% novel_contigs) %>%
    group_by(gene) %>%
    summarise(case_gene_count_wt = sum(num_reads_case, na.rm = TRUE), .groups = "drop")

#############################################################
# Calculate VAFs (TPM Based)
#############################################################
print("Estimating VAFs...")
# dropInfReps=TRUE must match the raw-count import above: without it, tximport may try to read
# Piscem/Oarfish Parquet bootstrap replicates, which requires the R package `arrow` in the container.
txi_tpm <- tximport(quant_file, type="oarfish", countsFromAbundance = "lengthScaledTPM",
                     txOut=TRUE, dropInfReps=TRUE)
qn <- data.frame(TPM = txi_tpm$abundance[, 1], transcript = clean_tx_names(rownames(txi_tpm$abundance)))

x <- inner_join(qn, tx2g, by = "transcript")
wt_tpm <- x[!x$transcript %in% novel_contigs, ] %>%
    group_by(gene) %>%
    summarise(WT = sum(TPM, na.rm = TRUE), .groups = "drop")

x <- left_join(x, wt_tpm, by = "gene")
x <- x %>% filter(transcript %in% cinfo$contig_id)

# Multi-gene variants (e.g. fusions):
# Rationale: average WT expression across partner genes so the denominator is not as large as
# sum(WT) (which would usually lower VAF). Use sum() instead if you want total WT pool across genes.
# - mean_WT_TPM matches the original script (mean WT across partner genes)
# - sum_WT_TPM is an alternative denominator (sum WT across partner genes)
wt_pools <- data.table(x)[,
                          list(
                              mean_WT_TPM = mean(WT, na.rm = TRUE),
                              sum_WT_TPM  = sum(WT, na.rm = TRUE)
                          ),
                          by = "transcript"]

x <- distinct(inner_join(x, wt_pools, by = "transcript"))

x$VAF <- x$TPM / (x$TPM + x$mean_WT_TPM)
x$VAF[is.nan(x$VAF)] <- 0

x$VAF_sum_WT_TPM <- x$TPM / (x$TPM + x$sum_WT_TPM)
x$VAF_sum_WT_TPM[is.nan(x$VAF_sum_WT_TPM)] <- 0

#############################################################
# Merge and Run Fisher's Exact Test
#############################################################
print("Running Statistical Enrichment Tests...")
x <- x %>%
    left_join(case_counts_raw, by = "transcript") %>%
    left_join(case_wt_gene_sum, by = "gene") %>%
    left_join(control_gene_total, by = "gene") %>%
    left_join(control_gene_wt, by = "gene") %>%
    left_join(control_counts_summary, by = "transcript") %>%
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# Control VAF and ratio (same as original script)
x$control_VAF <- x$control_contig_count / (x$control_contig_count + x$control_gene_count_wt)
x$control_VAF[is.nan(x$control_VAF) & x$control_contig_count == 0] <- 0
x$VAF_ratio <- x$VAF / x$control_VAF
x$VAF_ratio[is.infinite(x$VAF_ratio)] <- NA
x$VAF_ratio[is.nan(x$VAF_ratio)] <- NA

# Same 2x2 matrix for both tests: rows = [Variant, WT], cols = [Case, Control]
run_enrichment_tests <- function(r_case, wt_case, r_ctrl, wt_ctrl) {
    mat <- matrix(c(r_case, wt_case, r_ctrl, wt_ctrl), nrow = 2, byrow = TRUE)
    f <- fisher.test(mat, alternative = "greater")
    c(
        fisher_p_val = unname(f$p.value),
        odds_ratio = as.numeric(f$estimate)
    )
}

# Apply tests to each variant
if (nrow(x) > 0) {
    stats_results <- t(mapply(run_enrichment_tests,
                              x$num_reads_case, x$case_gene_count_wt,
                              x$control_contig_count, x$control_gene_count_wt))
    x$fisher_p_val <- stats_results[, "fisher_p_val"]
    x$odds_ratio <- stats_results[, "odds_ratio"]
} else {
    x$fisher_p_val <- numeric(0)
    x$odds_ratio <- numeric(0)
}

# Final formatting: original column order, then Fisher / case count columns
x <- x %>%
    rename(contig_id = transcript) %>%
    select(
        contig_id, gene, TPM, WT, mean_WT_TPM, sum_WT_TPM,
        control_gene_count_total, control_gene_count_wt,
        control_contig_count, control_VAF, VAF, VAF_ratio,
        VAF_sum_WT_TPM,
        fisher_p_val, odds_ratio,
        num_reads_case, case_gene_count_wt
    ) %>%
    arrange(fisher_p_val)

x$TPM <- signif(as.numeric(x$TPM), 2)
x$mean_WT_TPM <- signif(as.numeric(x$mean_WT_TPM), 2)
x$sum_WT_TPM <- signif(as.numeric(x$sum_WT_TPM), 2)
x$VAF <- signif(as.numeric(x$VAF), 2)
x$VAF_sum_WT_TPM <- signif(as.numeric(x$VAF_sum_WT_TPM), 2)
x$control_VAF <- signif(as.numeric(x$control_VAF), 2)
x$VAF_ratio <- signif(as.numeric(x$VAF_ratio), 2)
x$num_reads_case <- round(as.numeric(x$num_reads_case))
x$case_gene_count_wt <- round(as.numeric(x$case_gene_count_wt))
x$fisher_p_val <- signif(as.numeric(x$fisher_p_val), 2)
x$odds_ratio <- signif(as.numeric(x$odds_ratio), 2)

if (nrow(x) > 0) {
    write.table(x, file = outfile, row.names = FALSE, quote = FALSE, sep = "\t")
} else {
    cat("WARNING: No variants found for VAF estimation. Writing empty output file.\n")
    write.table(x, file = outfile, row.names = FALSE, quote = FALSE, sep = "\t")
}
print(paste("Results written to", outfile))