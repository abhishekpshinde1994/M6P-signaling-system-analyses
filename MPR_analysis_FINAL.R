# Enhanced MPR Sequence Similarity Analysis
# Complete working script - no verification steps
# Author: Abhishek Shinde

cat(paste(rep("=", 70), collapse=""), "\n")
cat("   MPR SEQUENCE SIMILARITY ANALYSIS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Load packages quietly
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(ape)
  library(seqinr)
  library(grid)
  library(gridExtra)
})

# Try to load plotly for interactive plots
plotly_available <- FALSE
tryCatch({
  suppressPackageStartupMessages(library(plotly))
  plotly_available <- TRUE
  cat("✓ Interactive plotting enabled (plotly loaded)\n")
}, error = function(e) {
  cat("! plotly not available - installing now...\n")
  install.packages("plotly", dependencies = TRUE, quiet = TRUE)
  tryCatch({
    suppressPackageStartupMessages(library(plotly))
    plotly_available <- TRUE
    cat("✓ plotly installed and loaded successfully\n")
  }, error = function(e2) {
    cat("! Could not install plotly - interactive plots will be skipped\n")
  })
})
cat("\n")

# 1. PARSE FASTA
cat("Parsing FASTA file...\n")
sequences <- read.fasta("MPR.txt", seqtype = "AA", as.string = TRUE)
sequences <- lapply(sequences, function(x) toupper(gsub("[^A-Z]", "", x[1])))
sequences <- sequences[sapply(sequences, nchar) > 20]
seq_names <- names(sequences)
cat("✓ Parsed", length(sequences), "sequences\n\n")

# 2. EXTRACT METADATA
extract_metadata <- function(seq_names) {
  metadata <- data.frame(
    ID = seq_names,
    Type = character(length(seq_names)),
    Organism = character(length(seq_names)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(seq_names)) {
    name <- seq_names[i]
    
    if (grepl("CIMPR", name)) {
      metadata$Type[i] <- "CI-MPR"
    } else if (grepl("CDMPR", name)) {
      metadata$Type[i] <- "CD-MPR"
    } else if (grepl("PMPR", name)) {
      metadata$Type[i] <- "P-MPR"
    } else {
      metadata$Type[i] <- "Unknown"
    }
    
    if (grepl("Trichomonas", name)) {
      metadata$Organism[i] <- "Trichomonas"
    } else if (grepl("Giardia", name)) {
      metadata$Organism[i] <- "Giardia"
    } else if (grepl("Entamoeba", name)) {
      metadata$Organism[i] <- "Entamoeba"
    } else if (grepl("Acanthamoeba", name)) {
      metadata$Organism[i] <- "Acanthamoeba"
    } else if (grepl("Dictyostelium", name)) {
      metadata$Organism[i] <- "Dictyostelium"
    } else if (grepl("Mastigamoeba", name)) {
      metadata$Organism[i] <- "Mastigamoeba"
    } else {
      metadata$Organism[i] <- "Unknown"
    }
  }
  return(metadata)
}

metadata <- extract_metadata(seq_names)

# 3. CALCULATE SIMILARITY MATRIX
cat("Calculating similarity matrix...\n")
calculate_identity <- function(seq1, seq2) {
  max_len <- max(nchar(seq1), nchar(seq2))
  min_len <- min(nchar(seq1), nchar(seq2))
  matches <- sum(strsplit(substr(seq1, 1, min_len), "")[[1]] == 
                   strsplit(substr(seq2, 1, min_len), "")[[1]])
  return(matches / max_len)
}

n <- length(sequences)
sim_matrix <- matrix(1, n, n)
rownames(sim_matrix) <- seq_names
colnames(sim_matrix) <- seq_names

for (i in 1:(n-1)) {
  if (i %% 10 == 0) cat("  Progress:", i, "/", n, "\n")
  for (j in (i+1):n) {
    sim <- calculate_identity(sequences[[i]], sequences[[j]])
    sim_matrix[i, j] <- sim
    sim_matrix[j, i] <- sim
  }
}
cat("✓ Calculation complete!\n\n")

# 4. BASIC STATISTICS
lengths <- sapply(sequences, nchar)
cat(paste(rep("=", 70), collapse=""), "\n")
cat("   DATASET SUMMARY\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Total sequences:", length(sequences), "\n")
cat("\nSequence lengths:\n")
cat("  Range:", min(lengths), "-", max(lengths), "aa\n")
cat("  Mean:", round(mean(lengths), 1), "aa\n")
cat("\nMPR Types:\n")
type_counts <- table(metadata$Type)
for (i in seq_along(type_counts)) {
  cat("  ", names(type_counts)[i], ":", type_counts[i], "\n")
}
cat("\nOrganisms:\n")
org_counts <- table(metadata$Organism)
for (i in seq_along(org_counts)) {
  cat("  ", names(org_counts)[i], ":", org_counts[i], "\n")
}

upper_tri <- sim_matrix[upper.tri(sim_matrix)]
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("   SIMILARITY STATISTICS\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Pairwise comparisons:", length(upper_tri), "\n")
cat(sprintf("Mean similarity: %.4f ± %.4f\n", mean(upper_tri), sd(upper_tri)))
cat(sprintf("Range: %.4f - %.4f\n\n", min(upper_tri), max(upper_tri)))

# 5. GROUP COMPARISONS
cat(paste(rep("=", 70), collapse=""), "\n")
cat("   GROUP-BASED SIMILARITY ANALYSIS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

within_type <- c()
between_type <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    type_i <- metadata$Type[metadata$ID == seq_names[i]]
    type_j <- metadata$Type[metadata$ID == seq_names[j]]
    if (type_i != "Unknown" && type_j != "Unknown") {
      if (type_i == type_j) {
        within_type <- c(within_type, sim_matrix[i, j])
      } else {
        between_type <- c(between_type, sim_matrix[i, j])
      }
    }
  }
}

if (length(within_type) > 0 && length(between_type) > 0) {
  cat("MPR Type Comparison:\n")
  cat(sprintf("  Within-type similarity: %.4f ± %.4f\n", mean(within_type), sd(within_type)))
  cat(sprintf("  Between-type similarity: %.4f ± %.4f\n", mean(between_type), sd(between_type)))
  test_result <- wilcox.test(within_type, between_type, alternative = "greater")
  cat(sprintf("  Mann-Whitney U test P-value: %.2e\n", test_result$p.value))
  if (test_result$p.value < 0.001) {
    cat("  *** HIGHLY SIGNIFICANT difference!\n")
  } else if (test_result$p.value < 0.05) {
    cat("  ** Significant difference\n")
  }
}

within_org <- c()
between_org <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    org_i <- metadata$Organism[metadata$ID == seq_names[i]]
    org_j <- metadata$Organism[metadata$ID == seq_names[j]]
    if (org_i != "Unknown" && org_j != "Unknown") {
      if (org_i == org_j) {
        within_org <- c(within_org, sim_matrix[i, j])
      } else {
        between_org <- c(between_org, sim_matrix[i, j])
      }
    }
  }
}

if (length(within_org) > 0 && length(between_org) > 0) {
  cat("\nOrganism Comparison:\n")
  cat(sprintf("  Within-organism similarity: %.4f ± %.4f\n", mean(within_org), sd(within_org)))
  cat(sprintf("  Between-organism similarity: %.4f ± %.4f\n", mean(between_org), sd(between_org)))
  test_result <- wilcox.test(within_org, between_org, alternative = "greater")
  cat(sprintf("  Mann-Whitney U test P-value: %.2e\n", test_result$p.value))
  if (test_result$p.value < 0.001) {
    cat("  *** HIGHLY SIGNIFICANT difference!\n")
  } else if (test_result$p.value < 0.05) {
    cat("  ** Significant difference\n")
  }
}
cat("\n")

# 6. EXTREME PAIRS
cat(paste(rep("=", 70), collapse=""), "\n")
cat("   TOP 10 MOST SIMILAR PAIRS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

upper_indices <- which(upper.tri(sim_matrix), arr.ind = TRUE)
similarities <- sim_matrix[upper_indices]
sorted_idx <- order(similarities, decreasing = TRUE)

similar_pairs <- data.frame()
for (k in 1:min(10, length(sorted_idx))) {
  idx <- sorted_idx[k]
  i <- upper_indices[idx, 1]
  j <- upper_indices[idx, 2]
  cat(sprintf("%d. Similarity: %.4f\n", k, sim_matrix[i, j]))
  cat("   ", substr(seq_names[i], 1, 60), "\n")
  cat("   (", metadata$Type[i], ", ", metadata$Organism[i], ")\n", sep = "")
  cat("   <->\n")
  cat("   ", substr(seq_names[j], 1, 60), "\n")
  cat("   (", metadata$Type[j], ", ", metadata$Organism[j], ")\n\n", sep = "")
  
  similar_pairs <- rbind(similar_pairs, data.frame(
    Seq1 = seq_names[i],
    Seq2 = seq_names[j],
    Similarity = round(sim_matrix[i, j], 4),
    Type1 = metadata$Type[i],
    Type2 = metadata$Type[j],
    Org1 = metadata$Organism[i],
    Org2 = metadata$Organism[j]
  ))
}

cat(paste(rep("=", 70), collapse=""), "\n")
cat("   TOP 10 MOST DIVERGENT PAIRS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

divergent_pairs <- data.frame()
for (k in 1:min(10, length(sorted_idx))) {
  idx <- sorted_idx[length(sorted_idx) - k + 1]
  i <- upper_indices[idx, 1]
  j <- upper_indices[idx, 2]
  cat(sprintf("%d. Similarity: %.4f\n", k, sim_matrix[i, j]))
  cat("   ", substr(seq_names[i], 1, 60), " (", metadata$Type[i], ")\n", sep = "")
  cat("   ", substr(seq_names[j], 1, 60), " (", metadata$Type[j], ")\n\n", sep = "")
  
  divergent_pairs <- rbind(divergent_pairs, data.frame(
    Seq1 = seq_names[i],
    Seq2 = seq_names[j],
    Similarity = round(sim_matrix[i, j], 4),
    Type1 = metadata$Type[i],
    Type2 = metadata$Type[j],
    Org1 = metadata$Organism[i],
    Org2 = metadata$Organism[j]
  ))
}

# 7. PER-SEQUENCE STATISTICS
seq_stats <- data.frame()
mean_sims <- sapply(1:n, function(i) mean(sim_matrix[i, ][sim_matrix[i, ] < 1]))
overall_mean <- mean(mean_sims)
overall_sd <- sd(mean_sims)

for (i in 1:n) {
  sims_no_self <- sim_matrix[i, ][sim_matrix[i, ] < 1]
  is_outlier <- abs(mean_sims[i] - overall_mean) > 2 * overall_sd
  
  seq_stats <- rbind(seq_stats, data.frame(
    ID = seq_names[i],
    Length = nchar(sequences[[i]]),
    Type = metadata$Type[i],
    Organism = metadata$Organism[i],
    Mean_Similarity = round(mean(sims_no_self), 4),
    Max_Similarity = round(max(sims_no_self), 4),
    Min_Similarity = round(min(sims_no_self), 4),
    SD_Similarity = round(sd(sims_no_self), 4),
    Is_Outlier = is_outlier
  ))
}

outliers <- seq_stats[seq_stats$Is_Outlier, ]
if (nrow(outliers) > 0) {
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("   OUTLIER SEQUENCES\n")
  cat(paste(rep("=", 70), collapse=""), "\n\n")
  cat("Found", nrow(outliers), "outlier sequence(s):\n\n")
  for (i in 1:nrow(outliers)) {
    cat("  ", substr(outliers$ID[i], 1, 70), "\n")
    cat("    Type:", outliers$Type[i], ", Organism:", outliers$Organism[i], "\n")
    cat("    Mean similarity:", outliers$Mean_Similarity[i], "\n\n")
  }
}

# 8. SAVE CSV FILES
cat("Saving CSV files...\n")
write.csv(seq_stats, "mpr_sequence_stats.csv", row.names = FALSE)
write.csv(sim_matrix, "mpr_similarity_matrix.csv")
write.csv(similar_pairs, "mpr_most_similar_pairs.csv", row.names = FALSE)
write.csv(divergent_pairs, "mpr_most_divergent_pairs.csv", row.names = FALSE)
write.csv(metadata, "mpr_metadata.csv", row.names = FALSE)
cat("  ✓ CSV files saved\n")

# 9. SPATIAL CLUSTERING ANALYSIS
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("   SPATIAL CLUSTERING ANALYSIS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# MDS (Multidimensional Scaling)
cat("Performing MDS analysis...\n")
dist_matrix <- as.dist(1 - sim_matrix)
mds_result <- cmdscale(dist_matrix, k = 2)
colnames(mds_result) <- c("MDS1", "MDS2")

# K-means clustering
optimal_k <- min(5, max(2, floor(n / 15)))
set.seed(42)
kmeans_result <- kmeans(mds_result, centers = optimal_k, nstart = 25)
clusters <- kmeans_result$cluster

cat("  Optimal clusters:", optimal_k, "\n")
cat("  Cluster sizes:\n")
for (i in 1:optimal_k) {
  cat("    Cluster", i, ":", sum(clusters == i), "sequences\n")
}

# PCA analysis
cat("\nPerforming PCA analysis...\n")
pca_result <- prcomp(sim_matrix, center = TRUE, scale. = TRUE)
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
cat(sprintf("  PC1 explains %.1f%% of variance\n", variance_explained[1]))
cat(sprintf("  PC2 explains %.1f%% of variance\n", variance_explained[2]))

# Save clustering results
cluster_data <- data.frame(
  ID = seq_names,
  Type = metadata$Type,
  Organism = metadata$Organism,
  Cluster = clusters,
  MDS1 = mds_result[, 1],
  MDS2 = mds_result[, 2],
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)
write.csv(cluster_data, "mpr_spatial_clusters.csv", row.names = FALSE)
cat("\n  ✓ Cluster data saved\n")

# 10. CREATE VISUALIZATIONS
cat("\nCreating visualizations...\n")

# Color schemes
type_colors <- c("CI-MPR" = "#E74C3C", "CD-MPR" = "#3498DB", "P-MPR" = "#2ECC71", "Unknown" = "#95A5A6")
cluster_colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C", "#E67E22", "#34495E")
row_colors <- type_colors[metadata$Type]
names(row_colors) <- seq_names

# === HEATMAP ===
ann_data <- data.frame(Type = metadata$Type, row.names = seq_names)
ann_colors <- list(Type = type_colors[names(type_colors) %in% unique(metadata$Type)])

png("mpr_similarity_heatmap.png", width = 3000, height = 3000, res = 300)
pheatmap(sim_matrix,
         color = colorRampPalette(c("#FFFFFF", "#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"))(100),
         annotation_row = ann_data,
         annotation_col = ann_data,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "MPR Sequence Similarity Matrix",
         border_color = NA)
dev.off()

pdf("mpr_similarity_heatmap.pdf", width = 10, height = 10)
pheatmap(sim_matrix,
         color = colorRampPalette(c("#FFFFFF", "#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"))(100),
         annotation_row = ann_data,
         annotation_col = ann_data,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "MPR Sequence Similarity Matrix",
         border_color = NA)
dev.off()
cat("  ✓ Heatmap saved (PNG + PDF)\n")

# === DISTRIBUTION PLOTS ===
png("mpr_distributions.png", width = 3000, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Panel 1: Overall distribution
hist(upper_tri, breaks = 50, col = "#3498DB", border = "black",
     main = "Overall Similarity Distribution", xlab = "Pairwise Similarity", ylab = "Count")
abline(v = mean(upper_tri), col = "red", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Mean = %.4f", mean(upper_tri)), 
       col = "red", lty = 2, lwd = 2, bty = "n")

# Panel 2: Within vs Between Type
if (length(within_type) > 0 && length(between_type) > 0) {
  boxplot(list(Within = within_type, Between = between_type),
          col = c("#2ECC71", "#E74C3C"),
          main = "Within vs Between Type",
          ylab = "Similarity Score")
}

# Panel 3: Length vs Similarity
plot(lengths, mean_sims, col = row_colors, pch = 16, cex = 1.2,
     main = "Length vs Mean Similarity",
     xlab = "Sequence Length (aa)", ylab = "Mean Similarity")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.7, bty = "n")

# Panel 4: Type distribution
type_table <- table(metadata$Type)
y_max <- max(type_table) * 1.2  # Y-axis 20% higher than tallest bar
barplot(type_table, col = type_colors[names(type_table)],
        main = "MPR Type Distribution", ylab = "Count", las = 2,
        ylim = c(0, y_max))

dev.off()

pdf("mpr_distributions.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(upper_tri, breaks = 50, col = "#3498DB", border = "black",
     main = "Overall Similarity Distribution", xlab = "Pairwise Similarity", ylab = "Count")
abline(v = mean(upper_tri), col = "red", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Mean = %.4f", mean(upper_tri)), 
       col = "red", lty = 2, lwd = 2, bty = "n")

if (length(within_type) > 0 && length(between_type) > 0) {
  boxplot(list(Within = within_type, Between = between_type),
          col = c("#2ECC71", "#E74C3C"),
          main = "Within vs Between Type",
          ylab = "Similarity Score")
}

plot(lengths, mean_sims, col = row_colors, pch = 16, cex = 1.2,
     main = "Length vs Mean Similarity",
     xlab = "Sequence Length (aa)", ylab = "Mean Similarity")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.7, bty = "n")

barplot(table(metadata$Type), col = type_colors[names(table(metadata$Type))],
        main = "MPR Type Distribution", ylab = "Count", las = 2,
        ylim = c(0, y_max))

dev.off()
cat("  ✓ Distribution plots saved (PNG + PDF)\n")

# === SPATIAL CLUSTERING MAPS ===
png("mpr_spatial_clustering.png", width = 3600, height = 2400, res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Panel 1: MDS by Type
plot(mds_result[, 1], mds_result[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = "MDS: Colored by MPR Type")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.8, bty = "n")

# Panel 2: MDS by Cluster
plot(mds_result[, 1], mds_result[, 2], 
     col = cluster_colors[clusters], pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = sprintf("MDS: K-means Clusters (k=%d)", optimal_k))
text(kmeans_result$centers[, 1], kmeans_result$centers[, 2], 
     labels = 1:optimal_k, cex = 1.5, font = 2, col = "black")

# Panel 3: MDS with labels for outliers
plot(mds_result[, 1], mds_result[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = "MDS: Outliers Labeled")
if (nrow(outliers) > 0) {
  outlier_idx <- which(seq_stats$Is_Outlier)
  text(mds_result[outlier_idx, 1], mds_result[outlier_idx, 2],
       labels = substr(seq_names[outlier_idx], 1, 15),
       pos = 3, cex = 0.6, col = "red")
}

# Panel 4: PCA by Type
plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
     main = "PCA: Colored by MPR Type")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.8, bty = "n")

# Label specific proteins
label_proteins <- c("CDMPR_Human", "CIMPR_Human", "TVMPR1", "TVMPR2", "TVMPR3", 
                    "TVMPR4", "TVMPR5", "TVMPR6", "TVMPR7")
for (pattern in label_proteins) {
  matches <- grep(pattern, seq_names, ignore.case = TRUE)
  if (length(matches) > 0) {
    for (idx in matches) {
      label_text <- if (grepl("Homo", pattern)) {
        ifelse(grepl("CDMPR", seq_names[idx]), "CDMPR_Human", "CIMPR_Human")
      } else {
        gsub(".*TVMPR([0-9]).*", "TVMPR\\1", seq_names[idx])
      }
      text(pca_result$x[idx, 1], pca_result$x[idx, 2], 
           labels = label_text, pos = 3, cex = 0.7, font = 2, col = "black")
    }
  }
}

# Panel 5: PCA by Cluster
plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = cluster_colors[clusters], pch = 16, cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
     main = sprintf("PCA: K-means Clusters (k=%d)", optimal_k))

# Label same specific proteins
for (pattern in label_proteins) {
  matches <- grep(pattern, seq_names, ignore.case = TRUE)
  if (length(matches) > 0) {
    for (idx in matches) {
      label_text <- if (grepl("Homo", pattern)) {
        ifelse(grepl("CDMPR", seq_names[idx]), "CDMPR_Human", "CIMPR_Human")
      } else {
        gsub(".*TVMPR([0-9]).*", "TVMPR\\1", seq_names[idx])
      }
      text(pca_result$x[idx, 1], pca_result$x[idx, 2], 
           labels = label_text, pos = 3, cex = 0.7, font = 2, col = "black")
    }
  }
}

# Panel 6: Cluster size distribution
cluster_table <- table(clusters)
cluster_y_max <- max(cluster_table) * 1.2
barplot(cluster_table, col = cluster_colors[1:optimal_k],
        main = "Cluster Size Distribution",
        xlab = "Cluster", ylab = "Number of Sequences",
        ylim = c(0, cluster_y_max))

dev.off()

pdf("mpr_spatial_clustering.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

plot(mds_result[, 1], mds_result[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = "MDS: Colored by MPR Type")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.8, bty = "n")

plot(mds_result[, 1], mds_result[, 2], 
     col = cluster_colors[clusters], pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = sprintf("MDS: K-means Clusters (k=%d)", optimal_k))
text(kmeans_result$centers[, 1], kmeans_result$centers[, 2], 
     labels = 1:optimal_k, cex = 1.5, font = 2, col = "black")

plot(mds_result[, 1], mds_result[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = "MDS Dimension 1", ylab = "MDS Dimension 2",
     main = "MDS: Outliers Labeled")
if (nrow(outliers) > 0) {
  outlier_idx <- which(seq_stats$Is_Outlier)
  text(mds_result[outlier_idx, 1], mds_result[outlier_idx, 2],
       labels = substr(seq_names[outlier_idx], 1, 15),
       pos = 3, cex = 0.6, col = "red")
}

plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = row_colors, pch = 16, cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
     main = "PCA: Colored by MPR Type")
legend("topright", legend = names(type_colors), 
       col = type_colors, pch = 16, cex = 0.8, bty = "n")

# Label specific proteins in PDF
label_proteins <- c("CDMPR_Human", "CIMPR_Human", "TVMPR1", "TVMPR2", "TVMPR3", 
                    "TVMPR4", "TVMPR5", "TVMPR6", "TVMPR7")
for (pattern in label_proteins) {
  matches <- grep(pattern, seq_names, ignore.case = TRUE)
  if (length(matches) > 0) {
    for (idx in matches) {
      label_text <- if (grepl("Homo", pattern)) {
        ifelse(grepl("CDMPR", seq_names[idx]), "CDMPR Human", "CIMPR Human")
      } else {
        gsub(".*TVMPR([0-9]).*", "TVMPR\\1", seq_names[idx])
      }
      text(pca_result$x[idx, 1], pca_result$x[idx, 2], 
           labels = label_text, pos = 3, cex = 0.8, font = 2, col = "black")
    }
  }
}

plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = cluster_colors[clusters], pch = 16, cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
     main = sprintf("PCA: K-means Clusters (k=%d)", optimal_k))

# Label same specific proteins in PDF
for (pattern in label_proteins) {
  matches <- grep(pattern, seq_names, ignore.case = TRUE)
  if (length(matches) > 0) {
    for (idx in matches) {
      label_text <- if (grepl("Homo", pattern)) {
        ifelse(grepl("CD-MPR", seq_names[idx]), "CD-MPR Human", "CI-MPR Human")
      } else {
        gsub(".*TVMPR([0-9]).*", "TVMPR\\1", seq_names[idx])
      }
      text(pca_result$x[idx, 1], pca_result$x[idx, 2], 
           labels = label_text, pos = 3, cex = 0.8, font = 2, col = "black")
    }
  }
}

barplot(cluster_table, col = cluster_colors[1:optimal_k],
        main = "Cluster Size Distribution",
        xlab = "Cluster", ylab = "Number of Sequences",
        ylim = c(0, cluster_y_max))

dev.off()
cat("  ✓ Spatial clustering maps saved (PNG + PDF)\n")

# === INTERACTIVE SPATIAL MAPS (HTML) ===
if (plotly_available) {
  cat("\nCreating interactive HTML plots...\n")
  
  # Prepare data for interactive plots
  plot_data <- data.frame(
    ID = seq_names,
    Short_Name = substr(seq_names, 1, 50),
    Type = metadata$Type,
    Organism = metadata$Organism,
    Length = lengths,
    Mean_Sim = round(mean_sims, 4),
    Cluster = clusters,
    MDS1 = mds_result[, 1],
    MDS2 = mds_result[, 2],
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Is_Outlier = seq_stats$Is_Outlier,
    stringsAsFactors = FALSE
  )
  
  # Color mappings
  type_color_map <- c("CI-MPR" = "#E74C3C", "CD-MPR" = "#3498DB", 
                      "P-MPR" = "#2ECC71", "Unknown" = "#95A5A6")
  cluster_color_map <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", 
                         "#9B59B6", "#1ABC9C", "#E67E22", "#34495E")
  
  # Interactive MDS plot colored by Type
  p_mds_type <- plot_ly(
    data = plot_data,
    x = ~MDS1,
    y = ~MDS2,
    color = ~Type,
    colors = type_color_map,
    text = ~paste0(
      "<b>", Short_Name, "</b><br>",
      "Type: ", Type, "<br>",
      "Organism: ", Organism, "<br>",
      "Length: ", Length, " aa<br>",
      "Mean Similarity: ", Mean_Sim, "<br>",
      "Cluster: ", Cluster,
      ifelse(Is_Outlier, "<br><b style='color:red;'>OUTLIER</b>", "")
    ),
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 10,
      line = list(color = "white", width = 1),
      opacity = 0.8
    )
  ) %>%
    layout(
      title = list(text = "<b>MDS Spatial Map - Colored by MPR Type</b><br><sub>Hover over points to see sequence details</sub>"),
      xaxis = list(title = "MDS Dimension 1", zeroline = FALSE),
      yaxis = list(title = "MDS Dimension 2", zeroline = FALSE),
      hovermode = "closest",
      plot_bgcolor = "#F5F5F5",
      paper_bgcolor = "white"
    )
  
  htmlwidgets::saveWidget(
    p_mds_type, 
    "mpr_interactive_mds_by_type.html",
    selfcontained = TRUE,
    title = "MPR MDS Analysis - By Type"
  )
  
  # Interactive MDS plot colored by Cluster
  p_mds_cluster <- plot_ly(
    data = plot_data,
    x = ~MDS1,
    y = ~MDS2,
    color = ~as.factor(Cluster),
    colors = cluster_color_map[1:optimal_k],
    text = ~paste0(
      "<b>", Short_Name, "</b><br>",
      "Type: ", Type, "<br>",
      "Organism: ", Organism, "<br>",
      "Length: ", Length, " aa<br>",
      "Mean Similarity: ", Mean_Sim, "<br>",
      "Cluster: ", Cluster,
      ifelse(Is_Outlier, "<br><b style='color:red;'>OUTLIER</b>", "")
    ),
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 10,
      line = list(color = "white", width = 1),
      opacity = 0.8
    )
  ) %>%
    layout(
      title = list(text = sprintf("<b>MDS Spatial Map - K-means Clusters (k=%d)</b><br><sub>Hover over points to see sequence details</sub>", optimal_k)),
      xaxis = list(title = "MDS Dimension 1", zeroline = FALSE),
      yaxis = list(title = "MDS Dimension 2", zeroline = FALSE),
      hovermode = "closest",
      plot_bgcolor = "#F5F5F5",
      paper_bgcolor = "white",
      showlegend = TRUE,
      legend = list(title = list(text = "Cluster"))
    )
  
  htmlwidgets::saveWidget(
    p_mds_cluster, 
    "mpr_interactive_mds_by_cluster.html",
    selfcontained = TRUE,
    title = "MPR MDS Analysis - By Cluster"
  )
  
  # Interactive PCA plot colored by Type
  p_pca_type <- plot_ly(
    data = plot_data,
    x = ~PC1,
    y = ~PC2,
    color = ~Type,
    colors = type_color_map,
    text = ~paste0(
      "<b>", Short_Name, "</b><br>",
      "Type: ", Type, "<br>",
      "Organism: ", Organism, "<br>",
      "Length: ", Length, " aa<br>",
      "Mean Similarity: ", Mean_Sim, "<br>",
      "Cluster: ", Cluster,
      ifelse(Is_Outlier, "<br><b style='color:red;'>OUTLIER</b>", "")
    ),
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 10,
      line = list(color = "white", width = 1),
      opacity = 0.8
    )
  ) %>%
    layout(
      title = list(text = "<b>PCA Spatial Map - Colored by MPR Type</b><br><sub>Hover over points to see sequence details</sub>"),
      xaxis = list(title = sprintf("PC1 (%.1f%% variance)", variance_explained[1]), zeroline = FALSE),
      yaxis = list(title = sprintf("PC2 (%.1f%% variance)", variance_explained[2]), zeroline = FALSE),
      hovermode = "closest",
      plot_bgcolor = "#F5F5F5",
      paper_bgcolor = "white"
    )
  
  htmlwidgets::saveWidget(
    p_pca_type, 
    "mpr_interactive_pca_by_type.html",
    selfcontained = TRUE,
    title = "MPR PCA Analysis - By Type"
  )
  
  # Interactive PCA plot colored by Cluster
  p_pca_cluster <- plot_ly(
    data = plot_data,
    x = ~PC1,
    y = ~PC2,
    color = ~as.factor(Cluster),
    colors = cluster_color_map[1:optimal_k],
    text = ~paste0(
      "<b>", Short_Name, "</b><br>",
      "Type: ", Type, "<br>",
      "Organism: ", Organism, "<br>",
      "Length: ", Length, " aa<br>",
      "Mean Similarity: ", Mean_Sim, "<br>",
      "Cluster: ", Cluster,
      ifelse(Is_Outlier, "<br><b style='color:red;'>OUTLIER</b>", "")
    ),
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 10,
      line = list(color = "white", width = 1),
      opacity = 0.8
    )
  ) %>%
    layout(
      title = list(text = sprintf("<b>PCA Spatial Map - K-means Clusters (k=%d)</b><br><sub>Hover over points to see sequence details</sub>", optimal_k)),
      xaxis = list(title = sprintf("PC1 (%.1f%% variance)", variance_explained[1]), zeroline = FALSE),
      yaxis = list(title = sprintf("PC2 (%.1f%% variance)", variance_explained[2]), zeroline = FALSE),
      hovermode = "closest",
      plot_bgcolor = "#F5F5F5",
      paper_bgcolor = "white",
      showlegend = TRUE,
      legend = list(title = list(text = "Cluster"))
    )
  
  htmlwidgets::saveWidget(
    p_pca_cluster, 
    "mpr_interactive_pca_by_cluster.html",
    selfcontained = TRUE,
    title = "MPR PCA Analysis - By Cluster"
  )
  
  # Combined interactive dashboard
  p_combined <- subplot(
    p_mds_type %>% layout(showlegend = TRUE),
    p_mds_cluster %>% layout(showlegend = TRUE),
    p_pca_type %>% layout(showlegend = TRUE),
    p_pca_cluster %>% layout(showlegend = TRUE),
    nrows = 2,
    shareX = FALSE,
    shareY = FALSE,
    titleX = TRUE,
    titleY = TRUE
  ) %>%
    layout(
      title = list(text = "<b>MPR Spatial Clustering Dashboard</b><br><sub>Interactive exploration of sequence relationships</sub>"),
      showlegend = FALSE
    )
  
  htmlwidgets::saveWidget(
    p_combined, 
    "mpr_interactive_dashboard.html",
    selfcontained = TRUE,
    title = "MPR Spatial Clustering Dashboard"
  )
  
  cat("  ✓ Interactive HTML plots saved:\n")
  cat("    - mpr_interactive_mds_by_type.html\n")
  cat("    - mpr_interactive_mds_by_cluster.html\n")
  cat("    - mpr_interactive_pca_by_type.html\n")
  cat("    - mpr_interactive_pca_by_cluster.html\n")
  cat("    - mpr_interactive_dashboard.html (all 4 plots combined)\n")
  cat("\n  ℹ️  Open HTML files in your web browser to explore!\n")
  cat("     (Hover over points to see sequence labels and details)\n")
} else {
  cat("\n! Skipping interactive plots (plotly not available)\n")
}

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("   ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("\nCSV Files:\n")
cat("  - mpr_sequence_stats.csv\n")
cat("  - mpr_similarity_matrix.csv\n")
cat("  - mpr_most_similar_pairs.csv\n")
cat("  - mpr_most_divergent_pairs.csv\n")
cat("  - mpr_metadata.csv\n")
cat("  - mpr_spatial_clusters.csv\n")
cat("\nStatic Visualizations (PNG + PDF):\n")
cat("  - mpr_similarity_heatmap.png/.pdf\n")
cat("  - mpr_distributions.png/.pdf\n")
cat("  - mpr_spatial_clustering.png/.pdf\n")
if (plotly_available) {
  cat("\nInteractive Visualizations (HTML):\n")
  cat("  - mpr_interactive_mds_by_type.html\n")
  cat("  - mpr_interactive_mds_by_cluster.html\n")
  cat("  - mpr_interactive_pca_by_type.html\n")
  cat("  - mpr_interactive_pca_by_cluster.html\n")
  cat("  - mpr_interactive_dashboard.html ⭐ (all 4 plots)\n")
  cat("\n  💡 Open HTML files in your browser to explore interactively!\n")
  cat("     Hover over points to see sequence labels and details.\n")
}
cat(paste(rep("=", 70), collapse=""), "\n")

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("   UNDERSTANDING THE SIMILARITY DISTRIBUTION PLOT\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("The 'Overall Similarity Distribution' plot shows:\n\n")

cat("📊 Y-AXIS = 'Count' (Frequency)\n")
cat("   This represents HOW MANY sequence pairs have a given similarity.\n\n")

cat("📊 X-AXIS = 'Pairwise Similarity' (0 to 1)\n")
cat("   This represents the similarity score between two sequences.\n\n")

cat("INTERPRETATION:\n")
cat("───────────────\n")
total_comparisons <- length(upper_tri)
cat(sprintf("• Total pairwise comparisons: %d\n", total_comparisons))
cat(sprintf("  (From %d sequences: %d × %d-1 / 2)\n\n", n, n, n))

cat("• Each bar in the histogram shows:\n")
cat("  'How many sequence pairs fall within this similarity range'\n\n")

cat("EXAMPLE from your data:\n")
# Get histogram data
h <- hist(upper_tri, breaks = 50, plot = FALSE)
peak_idx <- which.max(h$counts)
peak_count <- h$counts[peak_idx]
peak_range <- c(h$breaks[peak_idx], h$breaks[peak_idx + 1])

cat(sprintf("  Peak bar: ~%.0f sequence pairs have similarity between %.3f-%.3f\n", 
            peak_count, peak_range[1], peak_range[2]))
cat(sprintf("  This means: %.0f out of %d comparisons (%.1f%%) show this level\n",
            peak_count, total_comparisons, (peak_count/total_comparisons)*100))
cat("  of similarity.\n\n")

cat("WHY THIS MATTERS:\n")
cat("─────────────────\n")
cat("• Low similarity values (left side): Most of your sequences are here!\n")
cat(sprintf("  → Mean similarity = %.4f = only %.1f%% identity\n", 
            mean(upper_tri), mean(upper_tri)*100))
cat("  → This indicates ancient evolutionary divergence\n\n")

cat("• High similarity values (right side): Few sequence pairs\n")
cat(sprintf("  → Maximum similarity = %.4f = %.1f%% identity\n",
            max(upper_tri), max(upper_tri)*100))
cat("  → These are your most conserved sequences (likely orthologs)\n\n")

cat("• The RED dashed line shows the MEAN (average) similarity\n")
cat("  → It's your 'typical' similarity across all comparisons\n\n")

cat("WHAT THE SHAPE TELLS YOU:\n")
cat("─────────────────────────\n")
if (mean(upper_tri) < 0.3) {
  cat("• Your distribution is heavily LEFT-skewed\n")
  cat("  → Most pairs have LOW similarity (distant relationships)\n")
  cat("  → MPR protein family has undergone extensive divergence\n")
  cat("  → Different subfamilies (CI/CD/P) are very distinct\n")
} else {
  cat("• Your distribution shows moderate similarity\n")
  cat("  → Sequences are reasonably conserved\n")
  cat("  → Family members share recognizable similarity\n")
}

cat("\nBOTTOM LINE:\n")
cat("────────────\n")
cat(sprintf("Your %d sequences generated %d pairwise comparisons.\n", n, total_comparisons))
cat("The histogram shows that MOST of these comparisons (~", 
    round((sum(upper_tri < 0.1) / total_comparisons) * 100), "%)\n", sep="")
cat(sprintf("have very LOW similarity (<10%%), indicating these MPR sequences\n"))
cat(sprintf("are highly divergent from each other.\n\n"))

cat(paste(rep("=", 70), collapse=""), "\n\n")
