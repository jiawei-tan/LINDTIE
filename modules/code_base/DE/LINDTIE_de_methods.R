################################################################################
# Module      : LINDTIE_de_methods
# Description : run_edgeR function for differential expression analysis
# Copyright   : (c) Jia Wei Tan, Dec 2025
# License     : MIT
# Maintainer  : https://github.com/jiawei-tan
# Portability : POSIX
#
# Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
################################################################################

library(dplyr)
library(data.table)
library(edgeR)

# catchOarfish function (edgeR v4.8.0)
################################################################################
catchOarfish <- function(paths, verbose = TRUE)
#	Read transcriptwise counts and bootstrap samples from Oarfish output
#	Use bootstrap samples to estimate overdispersion of transcriptwise counts
#	Gordon Smyth
#	Created 4 July 2025. Last modified 4 July 2025.
{
	NSamples <- length(paths)

#	Use jsonlite and arrow packages for reading
	OK <- requireNamespace("jsonlite",quietly=TRUE)
	if(!OK) stop("jsonlite package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("readr",quietly=TRUE)
	if(!OK) stop("readr package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("arrow",quietly=TRUE)
	if(!OK) stop("arrow package required but is not installed (or can't be loaded)")

#	Initialize vector of inferential sample types
	ResampleType <- rep_len("bootstrap",NSamples)

#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",paths[j],", ", sep="")

#		File locations
		MetaFile <- paste0(paths[j],".meta_info.json")
		QuantFile <- paste0(paths[j],".quant")
		BootFile <- paste0(paths[j],".infreps.pq")
		if(!file.exists(QuantFile)) stop("quant file not found at specified path")

#		Meta information
		Meta <- jsonlite::fromJSON(MetaFile)
		NBoot <- Meta$num_bootstraps
		if(is.null(NBoot)) stop("Can't find number of bootstraps")
		if(verbose) cat(NBoot,"bootstraps\n")

#		Read counts
		if(j == 1L) {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdd",progress=FALSE))
			NTx <- nrow(Quant)
			Counts <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
			Counts[,1L] <- Quant$num_reads
			Ann <- data.frame(len=Quant$len)
			row.names(Ann) <- Quant$tname
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__d",progress=FALSE))
			Counts[,j] <- Quant$num_reads
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			Boot <- as.matrix(arrow::read_parquet(BootFile))
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}

#	Estimate overdispersion for each transcript
	i <- (DF > 0L)
	if(sum(i) > 0L) {
		OverDisp[i] <- OverDisp[i] / DF[i]
#		Apply a limited amount of moderation
		DFMedian <- median(DF[i])
		DFPrior <- 3
		OverDispPrior <- median(OverDisp[i]) / qf(0.5,df1=DFMedian,df2=DFPrior)
		if(OverDispPrior < 1) OverDispPrior <- 1
		OverDisp[i] <- (DFPrior * OverDispPrior + DF[i]*OverDisp[i]) / (DFPrior + DF[i])
		OverDisp <- pmax(OverDisp,1)
		OverDisp[!i] <- OverDispPrior
	} else {
		OverDisp[] <- NA_real_
		OverDispPrior <- NA_real_
	}

#	Prepare output
	dimnames(Counts) <- list(row.names(Ann),paths)
	Ann$Overdispersion <- OverDisp

	list(counts=Counts,annotation=Ann,overdispersion.prior=OverDispPrior,resample.type=ResampleType)
}

#########################################################################################

if (!exists("logfile")) {
  logfile <- "DE.log"
}
if (!exists("logcat")) {
  logcat <- function(...) {
    cat(..., file = logfile, append = TRUE)
    cat(...)
  }
}

run_edgeR <- function(case_name, oarfish_output_dir, outdir,
                      # cpm_cutoff = minCPM, qval = FDR,
                      # min_logfc = minLogFC) {
                      cpm_cutoff = NULL, qval = NULL,
                      min_logfc = NULL) {
  # Use global variables if provided, otherwise use defaults
  if (is.null(cpm_cutoff)) {
    if (exists("minCPM", envir = .GlobalEnv)) {
      cpm_cutoff <- get("minCPM", envir = .GlobalEnv)
    } else {
      cpm_cutoff <- 0.1
    }
  }
  if (is.null(qval)) {
    if (exists("FDR", envir = .GlobalEnv)) {
      qval <- get("FDR", envir = .GlobalEnv)
    } else {
      qval <- 0.05
    }
  }
  if (is.null(min_logfc)) {
    if (exists("minLogFC", envir = .GlobalEnv)) {
      min_logfc <- get("minLogFC", envir = .GlobalEnv)
    } else {
      min_logfc <- 2
    }
  }

  logcat("Starting LINDTIE_de_methods.R (run_edgeR) script...\n")
  logcat("Case name: ", case_name, "\n")
  logcat("Q-value (FDR)parameter: ", qval, "\n")
  logcat("CPM cutoff (minCPM) parameter: ", cpm_cutoff, "\n")
  logcat("Min logFC (minLogFC) parameter: ", min_logfc, "\n")

  quant_files <- list.files(oarfish_output_dir, ".quant$",
                            recursive = TRUE, full.names = TRUE)
  if (length(quant_files) == 0) {
    logcat("ERROR: No .quant files found in ", oarfish_output_dir, "\n")
    stop("No .quant files found in oarfish_output_dir")
  }
  logcat("Found ", length(quant_files), " quant files\n")

  quant_prefixes <- sub("\\.quant$", "", quant_files)
  logcat("Quant prefixes: ", paste(quant_prefixes, collapse = ", "), "\n")

  logcat("Calling catchOarfish...\n")
  catch <- catchOarfish(paths = quant_prefixes)
  logcat("catchOarfish completed successfully\n")
  logcat(
    "Catch object dimensions: ", nrow(catch$counts),
    " transcripts x ", ncol(catch$counts), " samples\n"
  )

  logcat("Bootstrap information:\n")
  for (i in 1:length(quant_prefixes)) {
    meta_file <- paste0(quant_prefixes[i], ".meta_info.json")
    if (file.exists(meta_file)) {
      meta <- jsonlite::fromJSON(meta_file)
      sample_name <- basename(quant_prefixes[i])
      n_boot <- meta$num_bootstraps
      logcat("  Sample ", sample_name, ": ", n_boot, " bootstraps\n")
    }
  }
  logcat("Resample type: ", paste(catch$resample.type, collapse = ", "), "\n")
  logcat("Overdispersion prior: ", catch$overdispersion.prior, "\n")
  logcat("Overdispersion summary:\n")
  logcat(capture.output(summary(catch$annotation$Overdispersion)), sep = "\n")

  divided_counts <- catch$counts/catch$annotation$Overdispersion
  logcat(
    "Divided counts dimensions: ", nrow(divided_counts),
    " transcripts x ", ncol(divided_counts), " samples\n"
  )
  logcat("Divided counts summary:\n")
  logcat(capture.output(summary(as.vector(divided_counts))), sep = "\n")
  
  # Extract sample names from full paths in column names
  sample_names <- basename(colnames(divided_counts))
  print(sample_names)
  logcat("Sample names: ", paste(sample_names, collapse = ", "), "\n")

  # Assign groups
  group <- rep("control", ncol(divided_counts))
  group[basename(colnames(divided_counts)) == case_name] <- "case"
  if (!"case" %in% group) {
    logcat("ERROR: No case sample found.\n")
    stop("No case sample found. Ensure case_name matches sample column.")
  }

  group <- factor(group, levels = c("control", "case"))
  logcat("Group assignments: ", paste(group, collapse = ", "), "\n")

  sample_info <- data.frame(SampleID = sample_names,
                            Group = group,
                            stringsAsFactors = FALSE)
  logcat("Sample info:\n")
  logcat(paste(capture.output(print(sample_info)), collapse = "\n"), "\n")

  # Create DGEList, filter, normalize
  dge <- DGEList(counts = divided_counts,
                 samples = sample_info,
                 group   = sample_info$Group)

  case_cpms <- cpm(dge)[, group == "case"] %>% as.numeric()

  # filter by CPM only
  keep_tx <- case_cpms > cpm_cutoff
  logcat(
    "Transcripts after filtering by CPM (", cpm_cutoff, "): ",
    sum(keep_tx), " retained, ", sum(!keep_tx), " filtered\n"
  )

  dge_filt  <- dge[keep_tx, , keep.lib.sizes = TRUE]

  # Check sample distribution
  n_controls <- sum(dge_filt$samples$Group == "control")
  n_cases <- sum(dge_filt$samples$Group == "case")
  logcat(
    "Sample distribution: ", n_controls, " controls, ",
    n_cases, " cases\n"
  )

  # Additional checks for data validity
  logcat("Checking data validity...\n")
  logcat(
    "Number of transcripts with zero counts in all samples: ",
    sum(rowSums(dge_filt$counts) == 0), "\n"
  )
  logcat(
    "Number of transcripts with non-zero counts: ",
    sum(rowSums(dge_filt$counts) > 0), "\n"
  )

  # Remove transcripts with zero counts in all samples
  keep_nonzero <- rowSums(dge_filt$counts) > 0
  if (sum(keep_nonzero) < nrow(dge_filt$counts)) {
    logcat(
      "Removing ", sum(!keep_nonzero),
      " transcripts with zero counts in all samples\n"
    )
    dge_filt <- dge_filt[keep_nonzero, , keep.lib.sizes = TRUE]
  }

  if (sum(keep_tx) == 0) {
    logcat("ERROR: All transcripts filtered out by cpm.\n")
    stop("All transcripts filtered out by cpm.")
  }

  # MDS plot
  # Open a PNG device
  png(filename = paste0(outdir, "/DE_MDS_plot.png"),
    width    = 800,
    height   = 600,
    res      = 150)
  plotMDS(dge_filt,
          col = as.numeric(dge_filt$samples$Group),
          main = "MDS Plot: Cases vs Controls",
          labels = dge_filt$samples$SampleID)
  legend("topright",  # Better position that won't overlap with plot
         legend = levels(dge_filt$samples$Group),
         col = 1:length(levels(dge_filt$samples$Group)),
         pch = 15)
  # Close the device
  dev.off()
  logcat("MDS plot written to DE_MDS_plot.png\n")

  # MD plot
  # Open a PNG device
  png(filename = paste0(outdir, "/DE_MD_plot.png"),
    width    = 800,
    height   = 600,
    res      = 150)
  plotMD(dge_filt,
          col = as.numeric(dge_filt$samples$Group),
          main = "MD Plot: Cases vs Controls") 
  # Close the device
  dev.off()
  logcat("MD plot written to DE_MD_plot.png\n")

  # Run DE analysis
  if (length(unique(group)) == 2 && table(group)["control"] > 1) {
    logcat("Using QLFit model for differential expression.\n")
    # Design matrix
    design <- model.matrix(~ Group, data = dge_filt$samples)

    # Dispersion estimation
    dge_disp <- estimateDisp(dge_filt, design = design, robust = TRUE)

    # Quasi-likelihood fitting and testing
    qlf_fit <- glmQLFit(dge_disp, design = design, robust = TRUE)
    qlf_res <- glmQLFTest(qlf_fit, coef = "Groupcase")

    # Extract results
    tryCatch({
      res_all <- topTags(qlf_res, n = Inf)$table
      logcat("Results extraction completed successfully\n")
    }, error = function(e) {
      logcat("ERROR in results extraction: ", e$message, "\n")
      stop("Failed to extract differential expression results")
    })

  } else {
    logcat("Only one control sample detected. ",
           "Using exact test with fixed dispersion.\n"
    )
    logcat("DGEList dimensions for exact test: ", nrow(dge$counts),
           " transcripts x ", ncol(dge$counts), " samples\n")
    et <- exactTest(dge, dispersion = 0.1)
    res_all <- as.data.frame(topTags(et, n = Inf))
    logcat("Exact test completed successfully\n")
  }

  logcat(
    "Found ", nrow(res_all), " differentially expressed transcripts\n"
  )

  # QL dispersion plot
  # Open a PNG device
  png(filename = paste0(outdir, "/DE_QLDisp_plot.png"),
    width    = 800,
    height   = 600,
    res      = 150)
  plotQLDisp(qlf_fit)
  dev.off()
  logcat("QL dispersion plot written to DE_QLDisp_plot.png\n")

  # Log FDR and logFC summaries BEFORE filtering
  logcat("FDR summary (before filtering):\n")
  logcat(capture.output(summary(res_all$FDR)), sep = "\n")
  logcat("logFC summary (before filtering):\n")
  logcat(capture.output(summary(res_all$logFC)), sep = "\n")

  # Apply filtering - Subset significant (based on FDR & logFC)
  res_sig <- subset(res_all, FDR < qval & logFC >= min_logfc)
  logcat("Filtering criteria: FDR < ", qval, " and logFC >= ", min_logfc, "\n")
  logcat("Found ", nrow(res_sig), " significantly differentially expressed transcripts (FDR < ", qval, ", logFC > ", min_logfc, ").\n")

  # remove CHESS
  res_sig <- res_sig[!grepl("CHS", rownames(res_sig)), ]
  logcat("Found ", nrow(res_sig), " significantly differentially expressed transcripts after removing reference\n")

if (nrow(res_sig) == 0) {
    logcat(
      "WARNING: No significant DE transcripts found after filtering ",
      "(FDR <", qval, ", logFC >", min_logfc, ").\n",
      "Proceeding with empty results.\n"
    )
  }

  # Add contig
  catch_anno <- catch$annotation
  catch_anno$transcript_id <- rownames(catch_anno)
  res_sig$transcript_id <- rownames(res_sig)
  logcat("Joining results with annotation...\n")
  res_sig <- left_join(res_sig, catch_anno, by = "transcript_id")
  res_sig <- res_sig[order(res_sig$PValue), ]

  # Add num_reads for the case sample to the results
  case_col <- which(sample_names == case_name)
  if (length(case_col) == 1) {
    # Get num_reads for the case sample
    num_reads_case <- catch$counts[, case_col]
    # Create a data.frame for joining - convert to integer
    num_reads_df <- data.frame(
      transcript_id = names(num_reads_case),
      num_reads_case = as.integer(num_reads_case),
      stringsAsFactors = FALSE
    )
    # Join by transcript_id
    res_sig <- left_join(res_sig, num_reads_df, by = "transcript_id")
  } else {
    logcat("WARNING: Could not find unique case sample for num_reads column.\n")
  }
  
  # Add total control reads across all control samples
  control_cols <- which(group == "control")
  if (length(control_cols) >= 1) {
    total_ctrl_reads <- rowSums(catch$counts[, control_cols, drop = FALSE])
    control_reads_df <- data.frame(
      transcript_id = rownames(catch$counts),
      total_num_reads_controls = as.integer(total_ctrl_reads),
      stringsAsFactors = FALSE
    )
    res_sig <- left_join(res_sig, control_reads_df, by = "transcript_id")
  } else {
    logcat("WARNING: No control samples found for total_num_reads_controls.\n")
  }
  
  # Round selected numeric columns to 2 significant figures
  cols_to_round <- c("logFC", "logCPM", "F", "PValue", "FDR")
  round_cols <- intersect(cols_to_round, names(res_sig))
  if (length(round_cols) > 0) {
    res_sig[round_cols] <- lapply(res_sig[round_cols], function(x) signif(as.numeric(x), 2))
  }
  logcat("Final results have", nrow(res_sig), " rows and ", ncol(res_sig), " columns\n")

  # Save filtered results
  out_file <- file.path(outdir, "DE_transcript_significant.txt")
  write.table(res_sig, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  logcat("Filtered significant DE results written to ", out_file, "\n")

  # Save full results (all transcripts)
  res_all$transcript_id <- rownames(res_all)
  res_all <- left_join(res_all, catch_anno, by = "transcript_id")
  res_all <- res_all[order(res_all$PValue), ]

  # Add num_reads for the case sample to the full results
  if (length(case_col) == 1) {
    res_all <- left_join(res_all, num_reads_df, by = "transcript_id")
  }

  # Add total control reads to full results
  if (length(control_cols) >= 1) {
    res_all <- left_join(res_all, control_reads_df, by = "transcript_id")
  }

  # Round selected numeric columns in full results to 2 significant figures
  round_cols_full <- intersect(cols_to_round, names(res_all))
  if (length(round_cols_full) > 0) {
    res_all[round_cols_full] <- lapply(res_all[round_cols_full], function(x) signif(as.numeric(x), 2))
  }

  out_file_full <- file.path(outdir, "DE_transcript_full_results.txt")
  write.table(res_all, file = out_file_full, sep = "\t", row.names = FALSE, quote = FALSE)
  logcat("Full DE results written to ", out_file_full, "\n")

  logcat("LINDTIE_de_methods.R (run_edgeR) script completed successfully.\n")

  return(res_sig)
}