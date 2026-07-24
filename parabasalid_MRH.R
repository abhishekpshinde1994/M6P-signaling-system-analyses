#!/usr/bin/env Rscript

################################################################################
# Human vs Parabasalid MRH Domain Similarity Analysis
################################################################################
# Author: Claude
# Purpose: Compare human MRH domains against T. vaginalis and Anaeramoeba
#          MRH domains to generate a 16x34 similarity matrix
################################################################################

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

required_packages <- c("ggplot2", "reshape2", "pheatmap", "RColorBrewer")

check_and_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("\n=== Checking required packages ===\n")
for (pkg in required_packages) {
  check_and_install(pkg)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
})

cat("✓ All required packages loaded successfully!\n")

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

#' Calculate pairwise sequence similarity (excluding gaps)
calculate_similarity_no_gaps <- function(seq1, seq2) {
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  if (length(s1) != length(s2)) {
    stop("Sequences must be of equal length (aligned)")
  }
  
  valid_pos <- which(s1 != "-" & s2 != "-")
  
  if (length(valid_pos) == 0) {
    return(0)
  }
  
  s1_valid <- s1[valid_pos]
  s2_valid <- s2[valid_pos]
  
  matches <- sum(s1_valid == s2_valid)
  similarity <- (matches / length(s1_valid)) * 100
  
  return(similarity)
}

#' Calculate pairwise sequence identity (including gaps)
calculate_identity_with_gaps <- function(seq1, seq2) {
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  if (length(s1) != length(s2)) {
    stop("Sequences must be of equal length (aligned)")
  }
  
  matches <- sum(s1 == s2)
  identity <- (matches / length(s1)) * 100
  
  return(identity)
}

#' Read FASTA file
read_fasta_base <- function(file_path) {
  cat(paste0("Reading FASTA file: ", file_path, "\n"))
  
  lines <- readLines(file_path)
  header_idx <- grep("^>", lines)
  
  sequences <- list()
  
  for (i in seq_along(header_idx)) {
    header <- lines[header_idx[i]]
    name <- sub("^>", "", header)
    
    if (i < length(header_idx)) {
      seq_lines <- lines[(header_idx[i] + 1):(header_idx[i + 1] - 1)]
    } else {
      seq_lines <- lines[(header_idx[i] + 1):length(lines)]
    }
    
    seq_lines <- seq_lines[seq_lines != ""]
    sequence <- paste(seq_lines, collapse = "")
    sequences[[name]] <- sequence
  }
  
  cat(paste0("✓ Read ", length(sequences), " sequences\n"))
  return(sequences)
}

#' Create rectangular similarity matrix
create_rectangular_similarity_matrix <- function(row_sequences, col_sequences, method = "no_gaps") {
  
  n_rows <- length(row_sequences)
  n_cols <- length(col_sequences)
  row_names <- names(row_sequences)
  col_names <- names(col_sequences)
  
  # Initialize matrix
  sim_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
  rownames(sim_matrix) <- row_names
  colnames(sim_matrix) <- col_names
  
  # Calculate pairwise similarities
  cat(paste0("Calculating ", n_rows, " x ", n_cols, " similarity matrix...\n"))
  pb <- txtProgressBar(min = 0, max = n_rows * n_cols, style = 3)
  counter <- 0
  
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      if (method == "no_gaps") {
        score <- calculate_similarity_no_gaps(row_sequences[[i]], col_sequences[[j]])
      } else {
        score <- calculate_identity_with_gaps(row_sequences[[i]], col_sequences[[j]])
      }
      sim_matrix[i, j] <- score
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  return(sim_matrix)
}

#' Clean sequence names for visualization
clean_name <- function(name) {
  # Extract key information for cleaner labels
  name <- gsub("_\\d+-\\d+.*$", "", name)  # Remove position info
  name <- gsub("_PDB:.*$", "", name)        # Remove PDB info
  return(name)
}

#' Get summary statistics
get_summary_stats <- function(matrix) {
  values <- as.vector(matrix)
  
  stats <- data.frame(
    Metric = c("Mean", "Median", "Min", "Max", "SD", "Q1", "Q3"),
    Value = c(
      mean(values),
      median(values),
      min(values),
      max(values),
      sd(values),
      quantile(values, 0.25),
      quantile(values, 0.75)
    )
  )
  
  return(stats)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

main <- function(fasta_file, output_dir = ".") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("HUMAN vs PARABASALID MRH DOMAIN SIMILARITY ANALYSIS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # ============================================================================
  # READ AND FILTER SEQUENCES
  # ============================================================================
  
  all_sequences <- read_fasta_base(fasta_file)
  
  # Filter sequences by category
  human_seqs <- list()
  parabasalid_seqs <- list()
  
  for (name in names(all_sequences)) {
    if (grepl("HUMAN", name, ignore.case = TRUE)) {
      human_seqs[[name]] <- all_sequences[[name]]
    } else if (grepl("TVMPR|Anaeramoeba", name, ignore.case = TRUE)) {
      parabasalid_seqs[[name]] <- all_sequences[[name]]
    }
  }
  
  cat("\n=== Sequence Categories ===\n")
  cat(paste0("Human MRH domains: ", length(human_seqs), "\n"))
  cat(paste0("Parabasalid MRH domains (TvMPR + Anaeramoeba): ", length(parabasalid_seqs), "\n"))
  
  if (length(human_seqs) == 0 || length(parabasalid_seqs) == 0) {
    stop("Error: Could not find sequences from both categories!")
  }
  
  # ============================================================================
  # CALCULATE SIMILARITY MATRICES
  # ============================================================================
  
  cat("\n=== Method 1: Similarity (ignoring gaps) ===\n")
  sim_matrix_nogaps <- create_rectangular_similarity_matrix(
    human_seqs, parabasalid_seqs, method = "no_gaps"
  )
  
  cat("\n=== Method 2: Identity (including gaps) ===\n")
  sim_matrix_gaps <- create_rectangular_similarity_matrix(
    human_seqs, parabasalid_seqs, method = "with_gaps"
  )
  
  # Clean names for visualization
  clean_row_names <- sapply(rownames(sim_matrix_nogaps), clean_name)
  clean_col_names <- sapply(colnames(sim_matrix_nogaps), clean_name)
  
  rownames(sim_matrix_nogaps) <- clean_row_names
  colnames(sim_matrix_nogaps) <- clean_col_names
  rownames(sim_matrix_gaps) <- clean_row_names
  colnames(sim_matrix_gaps) <- clean_col_names
  
  # ============================================================================
  # SAVE MATRICES
  # ============================================================================
  
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("SAVING SIMILARITY MATRICES\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  write.csv(sim_matrix_nogaps, 
            file.path(output_dir, "human_vs_parabasalid_similarity_no_gaps.csv"),
            row.names = TRUE)
  cat("✓ Saved: human_vs_parabasalid_similarity_no_gaps.csv\n")
  
  write.csv(sim_matrix_gaps, 
            file.path(output_dir, "human_vs_parabasalid_similarity_with_gaps.csv"),
            row.names = TRUE)
  cat("✓ Saved: human_vs_parabasalid_similarity_with_gaps.csv\n")
  
  save(sim_matrix_nogaps, sim_matrix_gaps, human_seqs, parabasalid_seqs,
       file = file.path(output_dir, "human_vs_parabasalid_matrices.RData"))
  cat("✓ Saved: human_vs_parabasalid_matrices.RData\n")
  
  # ============================================================================
  # GENERATE HEATMAPS
  # ============================================================================
  
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("GENERATING HEATMAPS\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # --- Heatmap 1: pheatmap style (no gaps) ---
  cat("Creating heatmap 1: Similarity (no gaps)...\n")
  
  pdf(file.path(output_dir, "heatmap_human_vs_parabasalid_no_gaps.pdf"), 
      width = 16, height = 10)
  
  pheatmap(sim_matrix_nogaps,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           display_numbers = TRUE,
           number_format = "%.1f",
           number_color = "black",
           fontsize_number = 6,
           fontsize = 9,
           fontsize_row = 9,
           fontsize_col = 8,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_method = "average",
           main = "Human vs Parabasalid MRH Domain Similarity (% - excluding gaps)",
           angle_col = 45,
           cellwidth = 15,
           cellheight = 18,
           border_color = "grey90",
           breaks = seq(0, 100, length.out = 101))
  
  dev.off()
  cat("✓ Saved: heatmap_human_vs_parabasalid_no_gaps.pdf\n")
  
  # --- Heatmap 2: pheatmap style (with gaps) ---
  cat("Creating heatmap 2: Identity (with gaps)...\n")
  
  pdf(file.path(output_dir, "heatmap_human_vs_parabasalid_with_gaps.pdf"), 
      width = 16, height = 10)
  
  pheatmap(sim_matrix_gaps,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           display_numbers = TRUE,
           number_format = "%.1f",
           number_color = "black",
           fontsize_number = 6,
           fontsize = 9,
           fontsize_row = 9,
           fontsize_col = 8,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_method = "average",
           main = "Human vs Parabasalid MRH Domain Identity (% - including gaps)",
           angle_col = 45,
           cellwidth = 15,
           cellheight = 18,
           border_color = "grey90",
           breaks = seq(0, 100, length.out = 101))
  
  dev.off()
  cat("✓ Saved: heatmap_human_vs_parabasalid_with_gaps.pdf\n")
  
  # --- Heatmap 3: No clustering (original order, no gaps) ---
  cat("Creating heatmap 3: No clustering (original order)...\n")
  
  pdf(file.path(output_dir, "heatmap_human_vs_parabasalid_no_clustering.pdf"), 
      width = 16, height = 10)
  
  pheatmap(sim_matrix_nogaps,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           display_numbers = TRUE,
           number_format = "%.1f",
           number_color = "black",
           fontsize_number = 6,
           fontsize = 9,
           fontsize_row = 9,
           fontsize_col = 8,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Human vs Parabasalid MRH Domain Similarity (original order)",
           angle_col = 45,
           cellwidth = 15,
           cellheight = 18,
           border_color = "grey90",
           breaks = seq(0, 100, length.out = 101))
  
  dev.off()
  cat("✓ Saved: heatmap_human_vs_parabasalid_no_clustering.pdf\n")
  
  # --- Heatmap 4: ggplot2 style ---
  cat("Creating heatmap 4: ggplot2 style...\n")
  
  melted_sim <- melt(sim_matrix_nogaps)
  colnames(melted_sim) <- c("Human_MRH", "Parabasalid_MRH", "Similarity")
  
  # Keep original order
  melted_sim$Human_MRH <- factor(melted_sim$Human_MRH, levels = rownames(sim_matrix_nogaps))
  melted_sim$Parabasalid_MRH <- factor(melted_sim$Parabasalid_MRH, levels = colnames(sim_matrix_nogaps))
  
  p1 <- ggplot(melted_sim, aes(x = Parabasalid_MRH, y = Human_MRH, fill = Similarity)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f", Similarity)), 
              size = 2, color = "black") +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")),
                         limits = c(0, 100),
                         name = "Similarity\n(%)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
    labs(x = "Parabasalid MRH Domains (TvMPR + Anaeramoeba)",
         y = "Human MRH Domains") +
    ggtitle("Human vs Parabasalid MRH Domain Similarity Matrix")
  
  ggsave(file.path(output_dir, "heatmap_human_vs_parabasalid_ggplot.pdf"), 
         plot = p1, width = 18, height = 10, units = "in")
  
  ggsave(file.path(output_dir, "heatmap_human_vs_parabasalid_ggplot.png"), 
         plot = p1, width = 18, height = 10, units = "in", dpi = 300)
  
  cat("✓ Saved: heatmap_human_vs_parabasalid_ggplot.pdf\n")
  cat("✓ Saved: heatmap_human_vs_parabasalid_ggplot.png\n")
  
  # --- Heatmap 5: Cluster columns only (rows in original order) ---
  cat("Creating heatmap 5: Cluster parabasalids only...\n")
  
  pdf(file.path(output_dir, "heatmap_human_vs_parabasalid_cluster_cols_only.pdf"), 
      width = 16, height = 10)
  
  pheatmap(sim_matrix_nogaps,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           display_numbers = TRUE,
           number_format = "%.1f",
           number_color = "black",
           fontsize_number = 6,
           fontsize = 9,
           fontsize_row = 9,
           fontsize_col = 8,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           clustering_method = "average",
           main = "Human MRH (rows) vs Parabasalid MRH (clustered columns)",
           angle_col = 45,
           cellwidth = 15,
           cellheight = 18,
           border_color = "grey90",
           breaks = seq(0, 100, length.out = 101))
  
  dev.off()
  cat("✓ Saved: heatmap_human_vs_parabasalid_cluster_cols_only.pdf\n")
  
  # ============================================================================
  # FIND TOP MATCHES
  # ============================================================================
  
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("FINDING TOP MATCHES\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # Find top 3 matches for each human domain
  top_matches <- data.frame()
  
  for (i in 1:nrow(sim_matrix_nogaps)) {
    human_name <- rownames(sim_matrix_nogaps)[i]
    scores <- sim_matrix_nogaps[i, ]
    top_idx <- order(scores, decreasing = TRUE)[1:3]
    
    for (rank in 1:3) {
      idx <- top_idx[rank]
      top_matches <- rbind(top_matches, data.frame(
        Human_MRH = human_name,
        Rank = rank,
        Parabasalid_MRH = colnames(sim_matrix_nogaps)[idx],
        Similarity = scores[idx]
      ))
    }
  }
  
  write.csv(top_matches, 
            file.path(output_dir, "top_matches_per_human_domain.csv"),
            row.names = FALSE)
  cat("✓ Saved: top_matches_per_human_domain.csv\n")
  
  # ============================================================================
  # GENERATE SUMMARY STATISTICS
  # ============================================================================
  
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("CALCULATING SUMMARY STATISTICS\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  stats_nogaps <- get_summary_stats(sim_matrix_nogaps)
  stats_gaps <- get_summary_stats(sim_matrix_gaps)
  
  summary_stats <- data.frame(
    Metric = stats_nogaps$Metric,
    Similarity_NoGaps = round(stats_nogaps$Value, 2),
    Identity_WithGaps = round(stats_gaps$Value, 2)
  )
  
  write.csv(summary_stats, 
            file.path(output_dir, "human_vs_parabasalid_summary_statistics.csv"),
            row.names = FALSE)
  cat("✓ Saved: human_vs_parabasalid_summary_statistics.csv\n")
  
  # Calculate per-human domain statistics
  human_stats <- data.frame()
  for (i in 1:nrow(sim_matrix_nogaps)) {
    human_name <- rownames(sim_matrix_nogaps)[i]
    scores <- sim_matrix_nogaps[i, ]
    
    human_stats <- rbind(human_stats, data.frame(
      Human_MRH = human_name,
      Mean_Similarity = mean(scores),
      Max_Similarity = max(scores),
      Min_Similarity = min(scores),
      SD = sd(scores)
    ))
  }
  
  write.csv(human_stats, 
            file.path(output_dir, "per_human_domain_statistics.csv"),
            row.names = FALSE)
  cat("✓ Saved: per_human_domain_statistics.csv\n")
  
  # ============================================================================
  # PRINT SUMMARY
  # ============================================================================
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("OVERALL SUMMARY STATISTICS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  print(summary_stats, row.names = FALSE)
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("TOP 5 HUMAN DOMAINS (by mean similarity)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  human_stats_sorted <- human_stats[order(human_stats$Mean_Similarity, decreasing = TRUE), ]
  print(head(human_stats_sorted, 5), row.names = FALSE)
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("OUTPUT FILES GENERATED\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  output_files <- c(
    "human_vs_parabasalid_similarity_no_gaps.csv",
    "human_vs_parabasalid_similarity_with_gaps.csv",
    "human_vs_parabasalid_matrices.RData",
    "heatmap_human_vs_parabasalid_no_gaps.pdf",
    "heatmap_human_vs_parabasalid_with_gaps.pdf",
    "heatmap_human_vs_parabasalid_no_clustering.pdf",
    "heatmap_human_vs_parabasalid_ggplot.pdf",
    "heatmap_human_vs_parabasalid_ggplot.png",
    "heatmap_human_vs_parabasalid_cluster_cols_only.pdf",
    "top_matches_per_human_domain.csv",
    "human_vs_parabasalid_summary_statistics.csv",
    "per_human_domain_statistics.csv"
  )
  
  for (f in output_files) {
    cat(paste0("  ✓ ", f, "\n"))
  }
  
  cat(paste0("\nAll files saved to: ", normalizePath(output_dir), "\n\n"))
  cat("Analysis complete! ✓\n\n")
  
  return(list(
    similarity_no_gaps = sim_matrix_nogaps,
    identity_with_gaps = sim_matrix_gaps,
    summary = summary_stats,
    human_stats = human_stats,
    top_matches = top_matches
  ))
}

# ==============================================================================
# PARSE COMMAND LINE ARGUMENTS AND RUN
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

fasta_file <- NULL
output_dir <- "."

if (length(args) >= 1) {
  fasta_file <- args[1]
}
if (length(args) >= 2) {
  output_dir <- args[2]
}

if (is.null(fasta_file)) {
  fasta_file <- "MRH_aligned.fasta"
  cat("\nNo input file specified. Usage:\n")
  cat("  Rscript human_vs_parabasalid_similarity.R <fasta_file> [output_dir]\n\n")
  cat("Using default: MRH_aligned.fasta\n\n")
}

if (!file.exists(fasta_file)) {
  stop(paste0("Error: File not found: ", fasta_file))
}

# Run analysis
results <- main(fasta_file, output_dir)

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Matrix dimensions: ", nrow(results$similarity_no_gaps), " rows (Human) x ", 
    ncol(results$similarity_no_gaps), " columns (Parabasalid)\n", sep = "")
cat(paste(rep("=", 80), collapse = ""), "\n\n")