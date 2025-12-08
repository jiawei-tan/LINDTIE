################################################################################
# Module      : LINDTIE_estimate_VAF
# Description : Estimates VAF from Oarfish quantification results
# Copyright   : (c) Jia Wei Tan, Dec 2025
# License     : MIT
# Maintainer  : https://github.com/jiawei-tan
# Portability : POSIX
#
# Usage:
# Rscript LINDTIE_estimate_VAF.R <transcript_counts_matrix_file> \
#   <quant_file> \
#   <contig_info_file> \
#   <tx_ref_fasta> \
#   <tx2gene_file> \
#   <outfile>
#
# Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
################################################################################

#############################################################
# Parse arguments
#############################################################

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0 | (length(args) < 6 & !"--help" %in% args)) {
    cat("
        Invalid arguments.
    ")
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Calculate rough VAF estimation from Oarfish quantification results

        Usage:
            Rscript LINDTIE_estimate_VAF.R \\
            <transcript_counts_matrix_file> \\
            <quant_file> \\
            <contig_info_file> \\
            <tx_ref_fasta> \\
            <tx2gene_file> \\
            <outfile>
    ")

    q(save="no")
}

library(tximport)
library(data.table)
library(dplyr)

options(stringsAsFactors = FALSE,
        error = function(){dump.frames("estimate_vaf", to.file = TRUE); q()})

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
txi <- tximport(quant_file, type="oarfish", countsFromAbundance = "lengthScaledTPM", txOut=TRUE, dropInfReps=TRUE)
transcript_counts_matrix <- fread(transcript_counts_matrix_file)
cinfo <- fread(contig_info_file)
tx2g <- fread(tx2gene_file, col.names = c("transcript", "gene"))
txs <- fread(tx_ref_fasta, header = FALSE, sep = "\n")

#############################################################
# Prepare data
#############################################################

print("Preparing data...")
# keep only variants of interest
cinfo <- cinfo[cinfo$variant_of_interest, ]

# get overlapping gene info for novel annotated contigs
c2g <- data.frame(contig_id = cinfo$contig_id, gene=cinfo$overlapping_genes)
split_genes <- sapply(c2g$gene, function(x){strsplit(x, "\\||:")})
c2g <- data.frame(transcript = rep(c2g$contig_id, sapply(split_genes, length)),
                  gene = unlist(split_genes))
c2g <- distinct(c2g)
c2g <- c2g[c2g$gene != "", ]

# rename transcript_id to transcript
colnames(transcript_counts_matrix)[colnames(transcript_counts_matrix) == "transcript_id"] <- "transcript"

# match non-novel contig ECs to genes and combine with novel contigs
# we do this so that we have a gene mapping for every contig and transcript
like_refseq <- transcript_counts_matrix$transcript %like% "hg38_ncbiRefSeq"
if(any(like_refseq)) {
    transcript_counts_matrix[like_refseq, "transcript"] <- sapply(transcript_counts_matrix$transcript[like_refseq],
                                                   function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]})
    # also fix tx2g to remove extra suffix from transcripts, e.g. NM_030642[_2]
    tx2g$transcript <- sapply(tx2g$transcript, function(x){paste(strsplit(x, "_")[[1]][1:2], collapse="_")})
    tx2g <- distinct(tx2g)
}

transcript2g <- right_join(transcript_counts_matrix, tx2g, by = "transcript")

tx2g <- transcript2g[, c("transcript", "gene")]
tx2g <- distinct(rbind(tx2g, c2g))

# get list of all reference transcripts
# we have to get these from the fasta as some wildtype transcripts
# may not be in the tx2gene reference because they lack gene names
txs <- txs[grep("^>", txs$V1), ]
txs <- sapply(txs$V1, function(x){strsplit(x, " ")[[1]][1]})
txs <- as.character(sapply(txs, gsub, pattern = ">", replacement = ""))
like_refseq <- txs %like% "hg38_ncbiRefSeq"
if(any(like_refseq)) {
    txs <- as.character(sapply(txs, function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]}))
}

# extract all novel contigs
# as in the DE step, get all contigs that have
# a transcript containing no reference transcripts
tx_transcript <- transcript_counts_matrix
tx_transcript$novel <- !tx_transcript$transcript %in% txs
novel_contig_transcript <- tx_transcript[, all(novel), by = "transcript"]
novel_contig_transcript <- unique(novel_contig_transcript$transcript[novel_contig_transcript$V1])
novel_contigs <- unique(tx_transcript$transcript[tx_transcript$transcript%in%novel_contig_transcript])

#############################################################
# Calculate VAFs
#############################################################

print("Estimating VAFs...")

# calculate the wildtype TPM by summing TPMs of all wildtype
# transcripts (anything that isn't a novel contig) per gene
qn <- data.frame(TPM = txi$abundance[, 1], transcript = rownames(txi$abundance))
like_refseq <- qn$transcript %like% "hg38_ncbiRefSeq"
if(any(like_refseq)) {
    qn[like_refseq, "transcript"] <- sapply(qn[like_refseq, "transcript"],
                                            function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]})
}
x <- inner_join(qn, tx2g, by = "transcript")
wt_count <- data.table(x[!x$transcript %in% novel_contigs,])
wt_count <- distinct(wt_count)[, list(WT = sum(TPM, na.rm = TRUE)), by = "gene"]
wt_count <- wt_count[wt_count$WT >= 0, ]

# now add the wildtype counts back to the quant table
# and extract only novel contigs
x <- left_join(x, wt_count, by = "gene")
x <- x[x$transcript %in% cinfo$contig_id, ]

# if contigs span multiple genes, we need to get the mean TPM
mean_tpm <- data.table(x)[, list(mean_WT_TPM = mean(WT, na.rm = TRUE)), by = "transcript"]
x <- distinct(inner_join(x, mean_tpm, by = c("transcript")))
x$VAF <- x$TPM / (x$TPM + x$mean_WT_TPM)
x$VAF[is.nan(x$VAF) & x$TPM == 0] <- 0

colnames(x)[2] <- "contig_id"
x <- x[, c("contig_id", "gene", "TPM", "WT", "mean_WT_TPM", "VAF")]
x$TPM <- signif(as.numeric(x$TPM), 2)
x$mean_WT_TPM <- signif(as.numeric(x$mean_WT_TPM), 2)
x$VAF <- signif(as.numeric(x$VAF), 2)
if (nrow(x) > 0) {
    write.table(x, file = outfile, row.names = FALSE, quote = FALSE, sep = "\t")
} else {
    stop("ERROR: no variants to output. Please double-check your reference files.")
}