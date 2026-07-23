# Protein Co-occurrence Analysis Script - 53 Species with PNG and PDF Outputs
# Statistical measures for protein co-occurrence patterns

# Load required libraries
library(readxl)
library(pheatmap)
library(ggplot2)
library(dplyr)

cat("=== PROTEIN CO-OCCURRENCE ANALYSIS - 53 SPECIES ===\n")
cat("Analyzing statistical co-occurrence patterns between proteins\n\n")

# =============================================================================
# DATA INPUT AND VALIDATION
# =============================================================================

# Read the data
cat("Reading Excel file...\n")
data <- read_excel("correlation.xlsx", sheet = "Sheet1")

cat("=== DATA VALIDATION ===\n")
cat("Raw data dimensions:", dim(data), "\n")

# Extract species information and protein data
species_info <- data.frame(
  Supergroup = data[[1]],  # Column A
  Species = data[[2]]      # Column B
)

# Show raw species info
cat("Raw species info (first 10 rows):\n")
print(head(species_info, 10))

# Remove any rows with missing species names
complete_rows <- !is.na(species_info$Species) & species_info$Species != "" & 
  !is.na(species_info$Supergroup) & species_info$Supergroup != ""

cat("Rows with valid species and supergroup names:", sum(complete_rows), "\n")

species_info <- species_info[complete_rows, ]
protein_data <- data[complete_rows, 3:15]  # Columns C-P
protein_names <- colnames(protein_data)

cat("Final species count:", nrow(species_info), "\n")
cat("Expected species count: 53\n")

# Verify we have exactly 53 species
if(nrow(species_info) != 53) {
  cat("WARNING: Expected 53 species but found", nrow(species_info), "\n")
  cat("Please check your data for missing or extra entries.\n\n")
  
  cat("Species found:\n")
  for(i in 1:nrow(species_info)) {
    cat(sprintf("%2d: %-20s (%s)\n", i, species_info$Species[i], species_info$Supergroup[i]))
  }
} else {
  cat("✓ Confirmed: 53 species found\n")
}

cat("\nProtein names (columns C-P):\n")
for(i in 1:length(protein_names)) {
  cat(sprintf("%2d: %s\n", i, protein_names[i]))
}

# Create matrices
paralog_matrix <- as.matrix(protein_data)
paralog_matrix[is.na(paralog_matrix)] <- 0
binary_matrix <- ifelse(paralog_matrix > 0, 1, 0)
rownames(binary_matrix) <- species_info$Species
colnames(binary_matrix) <- protein_names

n_species <- nrow(binary_matrix)
n_proteins <- ncol(binary_matrix)

cat("\nFinal matrix dimensions:", n_species, "species ×", n_proteins, "proteins\n")

# Display protein presence summary
cat("\n=== PROTEIN PRESENCE SUMMARY ===\n")
presence_counts <- colSums(binary_matrix)
for(i in 1:length(presence_counts)) {
  cat(sprintf("%-15s: %2d/%d species (%.1f%%)\n", 
              protein_names[i], presence_counts[i], n_species, 
              100 * presence_counts[i] / n_species))
}

# Verify totals
cat("\nTotal species in each calculation:", n_species, "\n")
cat("All calculations will use exactly", n_species, "species\n\n")

# =============================================================================
# ENHANCED CO-OCCURRENCE CALCULATION FUNCTION
# =============================================================================

calculate_cooccurrence_stats <- function(proteinA_data, proteinB_data, proteinA_name, proteinB_name, total_species) {
  
  # Verify input length
  if(length(proteinA_data) != total_species || length(proteinB_data) != total_species) {
    stop("Data length mismatch: expected", total_species, "species")
  }
  
  # 2x2 contingency table
  both_absent <- sum(proteinA_data == 0 & proteinB_data == 0)
  proteinB_only <- sum(proteinA_data == 0 & proteinB_data == 1)
  proteinA_only <- sum(proteinA_data == 1 & proteinB_data == 0)
  both_present <- sum(proteinA_data == 1 & proteinB_data == 1)
  
  # Verify totals
  total_check <- both_absent + proteinB_only + proteinA_only + both_present
  if(total_check != total_species) {
    cat("ERROR: Total mismatch for", proteinA_name, "vs", proteinB_name, 
        "- got", total_check, "expected", total_species, "\n")
  }
  
  # Basic percentages
  pct_both_present <- both_present / total_species * 100
  pct_both_absent <- both_absent / total_species * 100
  pct_only_one <- (proteinA_only + proteinB_only) / total_species * 100
  
  # 1. JACCARD INDEX (Primary co-occurrence measure)
  # J = |A ∩ B| / |A ∪ B| = Both_Present / (Both_Present + A_Only + B_Only)
  union_size <- both_present + proteinA_only + proteinB_only
  jaccard_index <- ifelse(union_size > 0, both_present / union_size, 0)
  
  # 2. DICE COEFFICIENT (Sørensen-Dice similarity)
  # D = 2|A ∩ B| / (|A| + |B|)
  total_A <- sum(proteinA_data)
  total_B <- sum(proteinB_data)
  dice_coeff <- ifelse((total_A + total_B) > 0, (2 * both_present) / (total_A + total_B), 0)
  
  # 3. CO-OCCURRENCE RATIO (Observed/Expected)
  expected_both <- (total_A * total_B) / total_species
  cooccurrence_ratio <- ifelse(expected_both > 0, both_present / expected_both, 0)
  
  # 4. PHI COEFFICIENT (Pearson correlation for binary data)
  numerator <- (both_present * both_absent) - (proteinA_only * proteinB_only)
  denominator <- sqrt((both_present + proteinA_only) * (both_present + proteinB_only) * 
                        (both_absent + proteinA_only) * (both_absent + proteinB_only))
  phi_coeff <- ifelse(denominator > 0, numerator / denominator, 0)
  
  # 5. YULE'S Q (Association measure)
  yules_numerator <- (both_present * both_absent) - (proteinA_only * proteinB_only)
  yules_denominator <- (both_present * both_absent) + (proteinA_only * proteinB_only)
  yules_q <- ifelse(yules_denominator > 0, yules_numerator / yules_denominator, 0)
  
  # 6. STATISTICAL SIGNIFICANCE TESTS
  # Create contingency table for tests
  cont_table <- matrix(c(both_absent, proteinB_only, proteinA_only, both_present), 
                       nrow = 2, byrow = TRUE)
  
  # Hypergeometric test (co-occurrence greater than chance)
  hypergeom_pvalue <- tryCatch({
    phyper(both_present - 1, total_A, total_species - total_A, total_B, lower.tail = FALSE)
  }, error = function(e) { 1.0 })
  
  # Fisher's exact test
  fisher_result <- tryCatch({
    fisher.test(cont_table, alternative = "greater")
  }, error = function(e) { list(p.value = 1.0, estimate = 1.0) })
  fisher_pvalue <- fisher_result$p.value
  fisher_odds_ratio <- fisher_result$estimate
  
  # Chi-square test
  chi_result <- tryCatch({
    chisq.test(cont_table, correct = FALSE)
  }, error = function(e) { list(statistic = 0, p.value = 1) })
  chi_square <- chi_result$statistic
  chi_pvalue <- chi_result$p.value
  cramers_v <- sqrt(chi_square / (total_species * min(1, 1)))  # For 2x2 table
  
  # 7. CO-OCCURRENCE STRENGTH CLASSIFICATION
  if(pct_both_present >= 70) {
    cooccurrence_strength <- "Very Strong"
  } else if(pct_both_present >= 50) {
    cooccurrence_strength <- "Strong"
  } else if(pct_both_present >= 30) {
    cooccurrence_strength <- "Moderate"
  } else if(pct_both_present >= 15) {
    cooccurrence_strength <- "Weak"
  } else {
    cooccurrence_strength <- "Very Weak"
  }
  
  # 8. RELATIONSHIP TYPE CLASSIFICATION
  if(jaccard_index >= 0.5 && pct_both_present >= 50) {
    relationship_type <- "Strong Co-occurrence"
  } else if(jaccard_index >= 0.3 && pct_both_present >= 30) {
    relationship_type <- "Moderate Co-occurrence"
  } else if(pct_only_one >= 50 && pct_both_present < 25) {
    relationship_type <- "Mutual Exclusivity"
  } else if(pct_both_absent >= 70) {
    relationship_type <- "Both Rare"
  } else {
    relationship_type <- "Independent"
  }
  
  return(list(
    proteinA = proteinA_name,
    proteinB = proteinB_name,
    total_species = total_species,
    both_absent = both_absent,
    proteinA_only = proteinA_only,
    proteinB_only = proteinB_only,
    both_present = both_present,
    union_size = union_size,
    pct_both_present = pct_both_present,
    pct_both_absent = pct_both_absent,
    pct_only_one = pct_only_one,
    jaccard_index = jaccard_index,
    dice_coefficient = dice_coeff,
    cooccurrence_ratio = cooccurrence_ratio,
    phi_coefficient = phi_coeff,
    yules_q = yules_q,
    hypergeom_pvalue = hypergeom_pvalue,
    fisher_pvalue = fisher_pvalue,
    fisher_odds_ratio = fisher_odds_ratio,
    chi_square = chi_square,
    chi_pvalue = chi_pvalue,
    cramers_v = cramers_v,
    cooccurrence_strength = cooccurrence_strength,
    relationship_type = relationship_type,
    significant_cooccurrence = (hypergeom_pvalue < 0.05 || fisher_pvalue < 0.05) && jaccard_index > 0.3
  ))
}

# =============================================================================
# ANALYZE ALL PROTEIN PAIRS WITH 53 SPECIES
# =============================================================================

cat("=== ANALYZING ALL PROTEIN PAIRS (53 SPECIES) ===\n")
cooccurrence_results <- data.frame()

for(i in 1:(n_proteins-1)) {
  for(j in (i+1):n_proteins) {
    
    proteinA_name <- protein_names[i]
    proteinB_name <- protein_names[j]
    
    proteinA_data <- binary_matrix[, i]
    proteinB_data <- binary_matrix[, j]
    
    # Calculate co-occurrence statistics with explicit species count
    stats <- calculate_cooccurrence_stats(proteinA_data, proteinB_data, 
                                          proteinA_name, proteinB_name, n_species)
    
    # Convert to data frame row
    result_row <- data.frame(
      ProteinA = stats$proteinA,
      ProteinB = stats$proteinB,
      Total_Species = stats$total_species,
      Both_Absent = stats$both_absent,
      ProteinA_Only = stats$proteinA_only,
      ProteinB_Only = stats$proteinB_only,
      Both_Present = stats$both_present,
      Union_Size = stats$union_size,
      Pct_Both_Present = round(stats$pct_both_present, 1),
      Pct_Both_Absent = round(stats$pct_both_absent, 1),
      Pct_Only_One = round(stats$pct_only_one, 1),
      Jaccard_Index = round(stats$jaccard_index, 4),
      Dice_Coefficient = round(stats$dice_coefficient, 4),
      Cooccurrence_Ratio = round(stats$cooccurrence_ratio, 3),
      Phi_Coefficient = round(stats$phi_coefficient, 4),
      Yules_Q = round(stats$yules_q, 4),
      Hypergeom_P = round(stats$hypergeom_pvalue, 6),
      Fisher_P = round(stats$fisher_pvalue, 6),
      Fisher_OR = round(stats$fisher_odds_ratio, 3),
      Chi_Square = round(stats$chi_square, 4),
      Chi_P = round(stats$chi_pvalue, 6),
      Cramers_V = round(stats$cramers_v, 4),
      Cooccurrence_Strength = stats$cooccurrence_strength,
      Relationship_Type = stats$relationship_type,
      Significant_Cooccurrence = stats$significant_cooccurrence,
      stringsAsFactors = FALSE
    )
    
    cooccurrence_results <- rbind(cooccurrence_results, result_row)
  }
}

cat("Completed analysis of", nrow(cooccurrence_results), "protein pairs using", n_species, "species\n\n")

# Verify all calculations used 53 species
total_checks <- unique(cooccurrence_results$Total_Species)
cat("Species totals in all calculations:", paste(total_checks, collapse = ", "), "\n")
if(length(total_checks) == 1 && total_checks[1] == 53) {
  cat("✓ Confirmed: All calculations used exactly 53 species\n\n")
} else {
  cat("⚠ WARNING: Inconsistent species totals found!\n\n")
}

# =============================================================================
# EXAMPLE VERIFICATION: CDMPR vs epsinR (53 SPECIES)
# =============================================================================

cat("=== EXAMPLE VERIFICATION: CDMPR vs epsinR (53 SPECIES) ===\n")

# Find CDMPR and epsinR
cdmpr_idx <- which(grepl("CDMPR|CD.MPR|CI.MPR", protein_names, ignore.case = TRUE))
epsinr_idx <- which(grepl("epsinR|epsin", protein_names, ignore.case = TRUE))

if(length(cdmpr_idx) > 0 && length(epsinr_idx) > 0) {
  cdmpr_name <- protein_names[cdmpr_idx[1]]
  epsinr_name <- protein_names[epsinr_idx[1]]
  
  cat("Found proteins:", cdmpr_name, "and", epsinr_name, "\n")
  
  # Find this pair in results
  pair_result <- cooccurrence_results[
    (cooccurrence_results$ProteinA == cdmpr_name & cooccurrence_results$ProteinB == epsinr_name) |
      (cooccurrence_results$ProteinA == epsinr_name & cooccurrence_results$ProteinB == cdmpr_name), ]
  
  if(nrow(pair_result) > 0) {
    cat("\n", cdmpr_name, "vs", epsinr_name, "Co-occurrence Analysis (53 species):\n")
    cat("Both present:", pair_result$Both_Present, "out of", pair_result$Total_Species, "species (", pair_result$Pct_Both_Present, "%)\n")
    cat(cdmpr_name, "only:", pair_result$ProteinA_Only, "species\n")
    cat(epsinr_name, "only:", pair_result$ProteinB_Only, "species\n")
    cat("Both absent:", pair_result$Both_Absent, "species\n")
    cat("Union size (species with at least one protein):", pair_result$Union_Size, "\n\n")
    
    cat("Jaccard Index Calculation:\n")
    cat("J =", pair_result$Both_Present, "/ (", pair_result$Both_Present, "+", 
        pair_result$ProteinA_Only, "+", pair_result$ProteinB_Only, ")\n")
    cat("J =", pair_result$Both_Present, "/", pair_result$Union_Size, "=", pair_result$Jaccard_Index, "\n\n")
    
    cat("Co-occurrence Statistics:\n")
    cat("Jaccard Index:", pair_result$Jaccard_Index, "(0-1 scale, higher = better co-occurrence)\n")
    cat("Dice Coefficient:", pair_result$Dice_Coefficient, "\n")
    cat("Co-occurrence Ratio:", pair_result$Cooccurrence_Ratio, "(>1 = higher than expected by chance)\n")
    cat("Phi Coefficient:", pair_result$Phi_Coefficient, "(-1 to +1, positive = co-occurrence)\n")
    cat("Fisher's Exact Test p-value:", pair_result$Fisher_P, "\n")
    cat("Hypergeometric Test p-value:", pair_result$Hypergeom_P, "\n")
    cat("Co-occurrence Strength:", pair_result$Cooccurrence_Strength, "\n")
    cat("Relationship Type:", pair_result$Relationship_Type, "\n")
    cat("Statistically Significant Co-occurrence:", pair_result$Significant_Cooccurrence, "\n\n")
  }
} else {
  cat("Could not find CDMPR or epsinR in the data\n")
}

# =============================================================================
# TOP CO-OCCURRING PROTEIN PAIRS
# =============================================================================

cat("=== TOP CO-OCCURRING PROTEIN PAIRS (53 SPECIES) ===\n")

# Sort by Jaccard Index
top_cooccurrence <- cooccurrence_results[order(cooccurrence_results$Jaccard_Index, decreasing = TRUE), ]
top_10_cooccurrence <- head(top_cooccurrence, 10)

cat("Top 10 protein pairs by co-occurrence (Jaccard Index):\n\n")
for(i in 1:nrow(top_10_cooccurrence)) {
  row <- top_10_cooccurrence[i, ]
  cat(sprintf("%d. %s vs %s:\n", i, row$ProteinA, row$ProteinB))
  cat(sprintf("   Both present: %d/%d species (%.1f%%)\n", row$Both_Present, row$Total_Species, row$Pct_Both_Present))
  cat(sprintf("   Union size: %d species, Jaccard Index: %.4f\n", row$Union_Size, row$Jaccard_Index))
  cat(sprintf("   Relationship: %s (%s)\n", row$Relationship_Type, row$Cooccurrence_Strength))
  cat(sprintf("   Fisher p=%.6f, Hypergeom p=%.6f\n\n", row$Fisher_P, row$Hypergeom_P))
}

# =============================================================================
# CREATE VISUALIZATION MATRICES (53 SPECIES)
# =============================================================================

cat("=== CREATING CO-OCCURRENCE HEATMAPS (53 SPECIES) ===\n")

# Initialize matrices
jaccard_matrix <- matrix(0, n_proteins, n_proteins)
cooccurrence_ratio_matrix <- matrix(1, n_proteins, n_proteins)
significance_matrix <- matrix(0, n_proteins, n_proteins)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- protein_names
rownames(cooccurrence_ratio_matrix) <- colnames(cooccurrence_ratio_matrix) <- protein_names
rownames(significance_matrix) <- colnames(significance_matrix) <- protein_names

# Fill matrices
for(i in 1:nrow(cooccurrence_results)) {
  row <- cooccurrence_results[i, ]
  protA_idx <- which(protein_names == row$ProteinA)
  protB_idx <- which(protein_names == row$ProteinB)
  
  # Symmetric matrices
  jaccard_matrix[protA_idx, protB_idx] <- jaccard_matrix[protB_idx, protA_idx] <- row$Jaccard_Index
  cooccurrence_ratio_matrix[protA_idx, protB_idx] <- cooccurrence_ratio_matrix[protB_idx, protA_idx] <- row$Cooccurrence_Ratio
  significance_matrix[protA_idx, protB_idx] <- significance_matrix[protB_idx, protA_idx] <- 
    ifelse(row$Significant_Cooccurrence, 1, 0)
}

# =============================================================================
# CREATE PNG AND PDF OUTPUTS (300 DPI)
# =============================================================================

cat("Creating PNG and PDF outputs at 300 DPI...\n")

# Function to create both PNG and PDF versions of each plot
create_dual_plots <- function(plot_data, plot_function, base_filename, title_text, ...) {
  
  # PNG version (300 DPI)
  png_filename <- paste0(base_filename, "_53species.png")
  png(png_filename, width = 16, height = 14, units = "in", res = 300)
  plot_function(plot_data, main = title_text, ...)
  dev.off()
  
  # PDF version
  pdf_filename <- paste0(base_filename, "_53species.pdf")
  pdf(pdf_filename, width = 16, height = 14)
  plot_function(plot_data, main = title_text, ...)
  dev.off()
  
  cat("Created:", png_filename, "and", pdf_filename, "\n")
}

# Jaccard Index Heatmap
cat("1. Creating Jaccard Index heatmaps...\n")
png("jaccard_cooccurrence_heatmap_53species.png", width = 16, height = 14, units = "in", res = 300)
pheatmap(jaccard_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "white",
         fontsize_number = 7,
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue", "navy"))(100),
         main = "Jaccard Co-occurrence Index (53 Species)\n(0 = no co-occurrence, 1 = perfect co-occurrence)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

pdf("jaccard_cooccurrence_heatmap_53species.pdf", width = 16, height = 14)
pheatmap(jaccard_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "white",
         fontsize_number = 7,
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue", "navy"))(100),
         main = "Jaccard Co-occurrence Index (53 Species)\n(0 = no co-occurrence, 1 = perfect co-occurrence)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

# Co-occurrence Ratio Heatmap
cat("2. Creating Co-occurrence Ratio heatmaps...\n")
png("cooccurrence_ratio_heatmap_53species.png", width = 16, height = 14, units = "in", res = 300)
pheatmap(cooccurrence_ratio_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "white",
         fontsize_number = 7,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Co-occurrence Ratio (53 Species)\n(<1 = less than expected, >1 = more than expected by chance)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

pdf("cooccurrence_ratio_heatmap_53species.pdf", width = 16, height = 14)
pheatmap(cooccurrence_ratio_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "white",
         fontsize_number = 7,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Co-occurrence Ratio (53 Species)\n(<1 = less than expected, >1 = more than expected by chance)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

# Statistical Significance Heatmap
cat("3. Creating Statistical Significance heatmaps...\n")
png("cooccurrence_significance_heatmap_53species.png", width = 16, height = 14, units = "in", res = 300)
pheatmap(significance_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         color = c("white", "darkgreen"),
         legend_breaks = c(0, 1),
         legend_labels = c("Not Significant", "Significant Co-occurrence"),
         main = "Statistical Significance of Co-occurrence (53 Species)\n(Green = Statistically Significant)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

pdf("cooccurrence_significance_heatmap_53species.pdf", width = 16, height = 14)
pheatmap(significance_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         color = c("white", "darkgreen"),
         legend_breaks = c(0, 1),
         legend_labels = c("Not Significant", "Significant Co-occurrence"),
         main = "Statistical Significance of Co-occurrence (53 Species)\n(Green = Statistically Significant)",
         fontsize = 9,
         fontsize_row = 7,
         fontsize_col = 7,
         border_color = "lightgray")
dev.off()

cat("All heatmaps created in both PNG (300 DPI) and PDF formats!\n")

# =============================================================================
# EXPORT RESULTS
# =============================================================================

cat("\n=== EXPORTING RESULTS (53 SPECIES) ===\n")

# Export detailed results
write.csv(cooccurrence_results, "detailed_cooccurrence_analysis_53species.csv", row.names = FALSE)
write.csv(jaccard_matrix, "jaccard_cooccurrence_matrix_53species.csv")
write.csv(cooccurrence_ratio_matrix, "cooccurrence_ratio_matrix_53species.csv")
write.csv(significance_matrix, "cooccurrence_significance_matrix_53species.csv")

# Export the binary matrix for verification
write.csv(binary_matrix, "binary_matrix_53species.csv")

cat("\n=== FILES CREATED ===\n")
cat("PNG FILES (300 DPI):\n")
cat("- jaccard_cooccurrence_heatmap_53species.png: Jaccard co-occurrence index heatmap\n")
cat("- cooccurrence_ratio_heatmap_53species.png: Co-occurrence ratio heatmap\n")
cat("- cooccurrence_significance_heatmap_53species.png: Statistical significance heatmap\n")

cat("\nPDF FILES:\n")
cat("- jaccard_cooccurrence_heatmap_53species.pdf: Jaccard co-occurrence index heatmap\n")
cat("- cooccurrence_ratio_heatmap_53species.pdf: Co-occurrence ratio heatmap\n")
cat("- cooccurrence_significance_heatmap_53species.pdf: Statistical significance heatmap\n")

cat("\nDATA FILES:\n")
cat("- detailed_cooccurrence_analysis_53species.csv: Complete co-occurrence statistics\n")
cat("- jaccard_cooccurrence_matrix_53species.csv: Jaccard index matrix\n")
cat("- cooccurrence_ratio_matrix_53species.csv: Co-occurrence ratio matrix\n")
cat("- cooccurrence_significance_matrix_53species.csv: Significance matrix\n")
cat("- binary_matrix_53species.csv: Binary presence/absence matrix for verification\n")

cat("\n=== FINAL SUMMARY (53 SPECIES) ===\n")
cat("Total species analyzed:", n_species, "\n")
cat("Total proteins analyzed:", n_proteins, "\n")
cat("Total protein pairs tested:", nrow(cooccurrence_results), "\n")
cat("Significant co-occurrences:", sum(cooccurrence_results$Significant_Cooccurrence), "\n")
cat("Mean Jaccard Index:", round(mean(cooccurrence_results$Jaccard_Index), 4), "\n")
cat("Max Jaccard Index:", round(max(cooccurrence_results$Jaccard_Index), 4), "\n")

# Relationship type summary
cat("\nRelationship Type Summary:\n")
relationship_summary <- table(cooccurrence_results$Relationship_Type)
for(i in 1:length(relationship_summary)) {
  cat(sprintf("%-25s: %d pairs\n", names(relationship_summary)[i], relationship_summary[i]))
}

cat("\n=== OUTPUT FORMAT SPECIFICATIONS ===\n")
cat("PNG FILES:\n")
cat("- Resolution: 300 DPI\n")
cat("- Dimensions: 16 × 14 inches\n")
cat("- Format: High-resolution PNG suitable for publication\n")

cat("\nPDF FILES:\n")
cat("- Dimensions: 16 × 14 inches\n")
cat("- Format: Vector-based PDF suitable for publication and printing\n")
cat("- Quality: Scalable without quality loss\n")

cat("\nAnalysis complete with exactly 53 species!\n")
cat("Both PNG (300 DPI) and PDF versions of all heatmaps have been created.\n")