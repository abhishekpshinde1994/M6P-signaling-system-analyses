# Proteomics IP Experiment Analysis - Perseus-style Pipeline
# Advanced Imputation, Student's T-test, and Fold Change Calculation
# Modified to include TVAG labels on volcano plots

# Load required libraries
library(readxl)      # For reading Excel files
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(tidyr)       # For data reshaping
library(plotly)      # For interactive plots
library(htmlwidgets) # For saving interactive plots
library(MASS)        # For multivariate normal distribution
library(ggrepel)     # For non-overlapping text labels in plots

# Check if ggrepel is installed and install if necessary
if (!require(ggrepel, quietly = TRUE)) {
  cat("Installing ggrepel package...\n")
  install.packages("ggrepel")
  library(ggrepel)
}

# ============================================================================
# CONFIGURATION - ADJUST THESE SETTINGS
# ============================================================================

# File path - CHANGE THIS TO YOUR FILE PATH
file_path <- "GNPT_DSP_IP.xlsx"

# Analysis parameters
imputation_method <- "perseus_downshift"  # Options: "perseus_downshift", "perseus_normal", "min_det", "knn", "remove_nan"
downshift_factor <- 1.8      # Width factor for downshifted normal distribution
downshift_mean_shift <- 0.3  # How much to shift mean down (as fraction of std)
min_detection_threshold <- 2 # Minimum replicates needed per group for reliable results
significance_threshold <- 0.05   # P-value threshold
fold_change_threshold <- 1.0     # Log2 fold change threshold
valid_values_threshold <- 2      # Minimum valid values needed in at least one group for imputation

# ============================================================================
# READ AND PREPARE DATA
# ============================================================================

cat("Reading proteomics data...\n")

# Read the Excel file
data <- read_excel(file_path, sheet = 1)

data <- data %>%
  mutate_all(~if_else(. == 0,
                      NA_character_,
                      .))

# Display data structure
cat("\nData structure:\n")
str(data)
cat("\nColumn names:\n")
print(colnames(data))
cat("\nFirst few rows:\n")
print(head(data, 3))

# Define column positions based on your structure:
# Columns 1-3: Protein annotations
# Columns 4-6: GNPT replicates (1834-GNP-DSP1/2/3)
# Columns 7-9: WT replicates (1834-WT-17-48-DSP1/2/3)

protein_genes <- data[[1]]        # T: PG.Genes
protein_descriptions <- data[[2]] # T: PG.ProteinDescriptions  
protein_uniprot <- data[[3]]      # T: PG.UniProtIds

GNPT_cols <- 4:6  # GNPT columns
wt_cols <- 7:9    # WT columns

# Extract data matrices
GNPT_original <- as.matrix(data[, GNPT_cols])
wt_original <- as.matrix(data[, wt_cols])

# Convert to numeric (in case Excel imported as character)
GNPT_original <- apply(GNPT_original, 2, as.numeric)
wt_original <- apply(wt_original, 2, as.numeric)

cat("\nOriginal data summary:\n")
cat("Total proteins:", nrow(data), "\n")
cat("GNPT replicates with NaN:", colSums(is.na(GNPT_original)), "\n")
cat("WT replicates with NaN:", colSums(is.na(wt_original)), "\n")
cat("Proteins with any GNPT NaN:", sum(rowSums(is.na(GNPT_original)) > 0), "\n")
cat("Proteins with any WT NaN:", sum(rowSums(is.na(wt_original)) > 0), "\n")
cat("Complete cases (no NaN):", sum(complete.cases(cbind(GNPT_original, wt_original))), "\n")

# ============================================================================
# PERSEUS-STYLE IMPUTATION FUNCTIONS
# ============================================================================

# Function to calculate distribution parameters for Perseus-style imputation
calculate_distribution_params <- function(data_matrix) {
  # Get all valid values
  valid_values <- as.vector(data_matrix)
  valid_values <- valid_values[!is.na(valid_values)]
  
  if (length(valid_values) < 3) {
    warning("Not enough valid values for distribution parameter estimation")
    return(list(mean = min(valid_values, na.rm = TRUE), sd = 1))
  }
  
  # Calculate parameters
  data_mean <- mean(valid_values)
  data_sd <- sd(valid_values)
  
  return(list(mean = data_mean, sd = data_sd))
}

# Perseus-style downshifted normal distribution imputation
perseus_downshift_imputation <- function(data_matrix, width_factor = 1.8, mean_shift = 0.3) {
  # Calculate distribution parameters from valid data
  params <- calculate_distribution_params(data_matrix)
  
  # Create downshifted distribution parameters
  # Mean is shifted down by mean_shift * standard deviation
  imputation_mean <- params$mean - mean_shift * params$sd
  # Width (standard deviation) is scaled by width_factor
  imputation_sd <- params$sd / width_factor
  
  cat("Perseus downshift imputation parameters:\n")
  cat("  Original mean:", round(params$mean, 3), "\n")
  cat("  Original SD:", round(params$sd, 3), "\n")
  cat("  Imputation mean:", round(imputation_mean, 3), "\n")
  cat("  Imputation SD:", round(imputation_sd, 3), "\n")
  
  # Apply imputation
  imputed_matrix <- data_matrix
  na_positions <- is.na(data_matrix)
  n_impute <- sum(na_positions)
  
  if (n_impute > 0) {
    # Generate random values from downshifted normal distribution
    imputed_values <- rnorm(n_impute, mean = imputation_mean, sd = imputation_sd)
    imputed_matrix[na_positions] <- imputed_values
  }
  
  return(imputed_matrix)
}

# Perseus-style normal distribution imputation (from valid values)
perseus_normal_imputation <- function(data_matrix) {
  # Calculate distribution parameters from valid data
  params <- calculate_distribution_params(data_matrix)
  
  cat("Perseus normal imputation parameters:\n")
  cat("  Mean:", round(params$mean, 3), "\n")
  cat("  SD:", round(params$sd, 3), "\n")
  
  # Apply imputation
  imputed_matrix <- data_matrix
  na_positions <- is.na(data_matrix)
  n_impute <- sum(na_positions)
  
  if (n_impute > 0) {
    # Generate random values from normal distribution of valid data
    imputed_values <- rnorm(n_impute, mean = params$mean, sd = params$sd)
    imputed_matrix[na_positions] <- imputed_values
  }
  
  return(imputed_matrix)
}

# Minimum detection value imputation (Perseus alterDSP)
min_detection_imputation <- function(data_matrix) {
  # Find minimum detected value
  min_detected <- min(data_matrix, na.rm = TRUE)
  
  # Use slightly lower value for imputation
  imputation_value <- min_detected - 0.1
  
  cat("Minimum detection imputation:\n")
  cat("  Minimum detected value:", round(min_detected, 3), "\n")
  cat("  Imputation value:", round(imputation_value, 3), "\n")
  
  # Apply imputation
  imputed_matrix <- data_matrix
  imputed_matrix[is.na(imputed_matrix)] <- imputation_value
  
  return(imputed_matrix)
}

# K-Nearest Neighbors imputation (simplified version)
knn_imputation <- function(data_matrix, k = 3) {
  cat("KNN imputation (k =", k, "):\n")
  
  imputed_matrix <- data_matrix
  
  # For each protein with missing values
  for (i in 1:nrow(data_matrix)) {
    if (any(is.na(data_matrix[i, ]))) {
      # Find k most similar proteins (by correlation)
      similarities <- numeric(nrow(data_matrix))
      
      for (j in 1:nrow(data_matrix)) {
        if (i != j) {
          # Calculate correlation using available values
          valid_i <- !is.na(data_matrix[i, ])
          valid_j <- !is.na(data_matrix[j, ])
          common_valid <- valid_i & valid_j
          
          if (sum(common_valid) >= 2) {
            similarities[j] <- cor(data_matrix[i, common_valid], 
                                   data_matrix[j, common_valid], 
                                   use = "complete.obs")
          } else {
            similarities[j] <- -1  # Low similarity for insufficient overlap
          }
        }
      }
      
      # Get k nearest neighbors
      similarities[i] <- -1  # Exclude self
      similarities[is.na(similarities)] <- -1
      nearest_indices <- order(similarities, decreasing = TRUE)[1:min(k, sum(similarities > -1))]
      
      # Impute missing values using mean of k nearest neighbors
      for (col in 1:ncol(data_matrix)) {
        if (is.na(data_matrix[i, col])) {
          neighbor_values <- data_matrix[nearest_indices, col]
          neighbor_values <- neighbor_values[!is.na(neighbor_values)]
          
          if (length(neighbor_values) > 0) {
            imputed_matrix[i, col] <- mean(neighbor_values)
          } else {
            # Fallback to minimum value if no neighbors have valid data
            imputed_matrix[i, col] <- min(data_matrix, na.rm = TRUE) - 0.1
          }
        }
      }
    }
  }
  
  return(imputed_matrix)
}

# Main imputation function with Perseus-style options
perform_perseus_imputation <- function(GNPT_matrix, wt_matrix, method = "perseus_downshift") {
  # Combine matrices to get global distribution parameters
  combined_matrix <- cbind(GNPT_matrix, wt_matrix)
  
  cat("\nPerforming Perseus-style imputation (method:", method, ")...\n")
  
  if (method == "perseus_downshift") {
    GNPT_imputed <- perseus_downshift_imputation(combined_matrix, 
                                                 width_factor = downshift_factor,
                                                 mean_shift = downshift_mean_shift)
    wt_imputed <- GNPT_imputed  # Use same imputation for consistency
    
    # Extract the respective columns
    GNPT_imputed <- GNPT_imputed[, 1:ncol(GNPT_matrix)]
    wt_imputed <- wt_imputed[, (ncol(GNPT_matrix)+1):ncol(combined_matrix)]
    
    # Apply original data where not missing
    GNPT_imputed[!is.na(GNPT_matrix)] <- GNPT_matrix[!is.na(GNPT_matrix)]
    wt_imputed[!is.na(wt_matrix)] <- wt_matrix[!is.na(wt_matrix)]
    
  } else if (method == "perseus_normal") {
    GNPT_imputed <- perseus_normal_imputation(combined_matrix)
    wt_imputed <- GNPT_imputed
    
    # Extract the respective columns
    GNPT_imputed <- GNPT_imputed[, 1:ncol(GNPT_matrix)]
    wt_imputed <- wt_imputed[, (ncol(GNPT_matrix)+1):ncol(combined_matrix)]
    
    # Apply original data where not missing
    GNPT_imputed[!is.na(GNPT_matrix)] <- GNPT_matrix[!is.na(GNPT_matrix)]
    wt_imputed[!is.na(wt_matrix)] <- wt_matrix[!is.na(wt_matrix)]
    
  } else if (method == "min_det") {
    GNPT_imputed <- min_detection_imputation(GNPT_matrix)
    wt_imputed <- min_detection_imputation(wt_matrix)
    
  } else if (method == "knn") {
    # For KNN, impute the combined matrix
    combined_imputed <- knn_imputation(combined_matrix, k = 3)
    GNPT_imputed <- combined_imputed[, 1:ncol(GNPT_matrix)]
    wt_imputed <- combined_imputed[, (ncol(GNPT_matrix)+1):ncol(combined_matrix)]
    
  } else if (method == "remove_nan") {
    cat("No imputation - using original data with NaN handling in statistics\n")
    return(list(GNPT = GNPT_matrix, wt = wt_matrix))
  }
  
  return(list(GNPT = GNPT_imputed, wt = wt_imputed))
}

# ============================================================================
# PERFORM IMPUTATION WITH PERSEUS METHODS
# ============================================================================

cat("\nPerforming Perseus-style imputation...\n")

# Set seed for reproducibility
set.seed(123)

# Filter proteins that have sufficient valid values for meaningful analysis
sufficient_data_mask <- (rowSums(!is.na(GNPT_original)) >= 1 | rowSums(!is.na(wt_original)) >= 1) &
  (rowSums(!is.na(GNPT_original)) + rowSums(!is.na(wt_original)) >= valid_values_threshold)

cat("Proteins with sufficient data for analysis:", sum(sufficient_data_mask), "out of", nrow(data), "\n")

# Filter data to proteins with sufficient valid values
GNPT_filtered <- GNPT_original[sufficient_data_mask, ]
wt_filtered <- wt_original[sufficient_data_mask, ]
protein_genes_filtered <- protein_genes[sufficient_data_mask]
protein_descriptions_filtered <- protein_descriptions[sufficient_data_mask]
protein_uniprot_filtered <- protein_uniprot[sufficient_data_mask]

# Perform imputation
imputation_result <- perform_perseus_imputation(GNPT_filtered, wt_filtered, method = imputation_method)
GNPT_imputed <- imputation_result$GNPT
wt_imputed <- imputation_result$wt

# Track imputation usage
GNPT_imputation_mask <- is.na(GNPT_filtered)
wt_imputation_mask <- is.na(wt_filtered)

cat("Imputation completed.\n")
cat("GNPT values imputed:", sum(GNPT_imputation_mask), "\n")
cat("WT values imputed:", sum(wt_imputation_mask), "\n")

# ============================================================================
# STATISTICAL ANALYSIS FUNCTIONS (Perseus-style)
# ============================================================================

# Function to calculate statistics for each protein (Perseus-style)
calculate_perseus_stats <- function(GNPT_vals, wt_vals, 
                                    GNPT_orig, wt_orig,
                                    min_detection = 2) {
  
  # Count valid original values
  valid_GNPT_orig <- sum(!is.na(GNPT_orig))
  valid_wt_orig <- sum(!is.na(wt_orig))
  
  # Check if imputation was used
  imputation_used <- any(is.na(GNPT_orig)) || any(is.na(wt_orig))
  
  # Calculate means
  GNPT_mean <- mean(GNPT_vals, na.rm = TRUE)
  wt_mean <- mean(wt_vals, na.rm = TRUE)
  
  # Calculate difference and fold change
  difference <- GNPT_mean - wt_mean  # This is log2 fold change since data is in log2
  log2_fc <- difference
  
  # Initialize statistical results
  t_stat <- NA
  p_value <- NA
  GNPT_sd <- NA
  wt_sd <- NA
  
  # Determine which data to use for statistics (Perseus philosophy)
  # Use imputed data if we don't have enough original values
  use_GNPT <- if (valid_GNPT_orig >= min_detection) GNPT_orig[!is.na(GNPT_orig)] else GNPT_vals
  use_wt <- if (valid_wt_orig >= min_detection) wt_orig[!is.na(wt_orig)] else wt_vals
  
  # Remove any remaining NAs
  use_GNPT <- use_GNPT[!is.na(use_GNPT)]
  use_wt <- use_wt[!is.na(use_wt)]
  
  # Perform t-test if we have enough data
  tryCatch({
    if (length(use_GNPT) >= 2 && length(use_wt) >= 2) {
      t_result <- t.test(use_GNPT, use_wt, var.equal = FALSE)
      t_stat <- as.numeric(t_result$statistic)
      p_value <- t_result$p.value
      GNPT_sd <- sd(use_GNPT)
      wt_sd <- sd(use_wt)
    } else if (length(use_GNPT) >= 1 && length(use_wt) >= 1) {
      # For single values, set p-value to 1 (no significance)
      GNPT_sd <- if(length(use_GNPT) > 1) sd(use_GNPT) else 0
      wt_sd <- if(length(use_wt) > 1) sd(use_wt) else 0
      p_value <- 1.0
      t_stat <- 0
    }
  }, error = function(e) {
    cat("Error in t-test calculation:", e$message, "\n")
  })
  
  # Determine data quality based on Perseus criteria
  data_quality <- case_when(
    valid_GNPT_orig >= min_detection && valid_wt_orig >= min_detection ~ "High",
    (valid_GNPT_orig + valid_wt_orig) >= min_detection ~ "Medium",
    TRUE ~ "Low"
  )
  
  # Determine imputation category
  imputation_category <- case_when(
    !imputation_used ~ "No_Imputation",
    valid_GNPT_orig == 0 || valid_wt_orig == 0 ~ "Full_Group_Imputed",
    TRUE ~ "Partial_Imputation"
  )
  
  return(list(
    GNPT_mean = GNPT_mean,
    wt_mean = wt_mean,
    difference = difference,
    log2_fc = log2_fc,
    t_statistic = t_stat,
    p_value = p_value,
    GNPT_sd = GNPT_sd,
    wt_sd = wt_sd,
    valid_GNPT_original = valid_GNPT_orig,
    valid_wt_original = valid_wt_orig,
    imputation_used = imputation_used,
    imputation_category = imputation_category,
    data_quality = data_quality
  ))
}

# ============================================================================
# PERFORM STATISTICAL ANALYSIS
# ============================================================================

cat("\nPerforming Perseus-style statistical analysis...\n")

# Initialize results dataframe
results <- data.frame(
  Protein_Gene = protein_genes_filtered,
  Protein_Description = protein_descriptions_filtered,
  UniProt_ID = protein_uniprot_filtered,
  GNPT_Mean = NA,
  WT_Mean = NA,
  Difference = NA,           # Mean difference (GNPT - WT)
  Log2_FC = NA,             # Log2 fold change
  T_Statistic = NA,
  P_Value = NA,
  Adj_P_Value = NA,         # FDR corrected
  GNPT_SD = NA,
  WT_SD = NA,
  Valid_GNPT_Original = NA,
  Valid_WT_Original = NA,
  Imputation_Used = FALSE,
  Imputation_Category = NA,
  Data_Quality = NA,
  stringsAsFactors = FALSE
)

# Progress bar
pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)

# Analyze each protein
for (i in 1:nrow(results)) {
  # Get values for current protein
  GNPT_vals <- GNPT_imputed[i, ]
  wt_vals <- wt_imputed[i, ]
  GNPT_orig <- GNPT_filtered[i, ]
  wt_orig <- wt_filtered[i, ]
  
  # Calculate statistics
  stats <- calculate_perseus_stats(GNPT_vals, wt_vals, GNPT_orig, wt_orig, 
                                   min_detection_threshold)
  
  # Store results
  results$GNPT_Mean[i] <- stats$GNPT_mean
  results$WT_Mean[i] <- stats$wt_mean
  results$Difference[i] <- stats$difference
  results$Log2_FC[i] <- stats$log2_fc
  results$T_Statistic[i] <- stats$t_statistic
  results$P_Value[i] <- stats$p_value
  results$GNPT_SD[i] <- stats$GNPT_sd
  results$WT_SD[i] <- stats$wt_sd
  results$Valid_GNPT_Original[i] <- stats$valid_GNPT_original
  results$Valid_WT_Original[i] <- stats$valid_wt_original
  results$Imputation_Used[i] <- stats$imputation_used
  results$Imputation_Category[i] <- stats$imputation_category
  results$Data_Quality[i] <- stats$data_quality
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Apply multiple testing correction (Benjamini-Hochberg FDR)
results$Adj_P_Value <- p.adjust(results$P_Value, method = "fdr")

cat("\nPereus-style statistical analysis completed.\n")

# ============================================================================
# FILTER AND CATEGORIZE RESULTS
# ============================================================================

cat("\nFiltering and categorizing results...\n")

# Complete results (with valid p-values)
complete_results <- results[!is.na(results$P_Value), ]

# High quality results
high_quality_results <- complete_results[complete_results$Data_Quality == "High", ]

# Medium quality results
medium_quality_results <- complete_results[complete_results$Data_Quality == "Medium", ]

# Significant results - updated to match volcano plot logic
significant_results <- complete_results[
  !is.na(complete_results$P_Value) & 
    complete_results$P_Value < significance_threshold &
    abs(complete_results$Log2_FC) > fold_change_threshold,
]

# Enriched proteins - updated to match volcano plot quadrants
enriched_GNPT <- complete_results[
  !is.na(complete_results$P_Value) & 
    complete_results$P_Value < significance_threshold &
    complete_results$Log2_FC > fold_change_threshold,
]

enriched_wt <- complete_results[
  !is.na(complete_results$P_Value) & 
    complete_results$P_Value < significance_threshold &
    complete_results$Log2_FC < -fold_change_threshold,
]

# Proteins by imputation category
no_imputation <- complete_results[complete_results$Imputation_Category == "No_Imputation", ]
partial_imputation <- complete_results[complete_results$Imputation_Category == "Partial_Imputation", ]
full_group_imputed <- complete_results[complete_results$Imputation_Category == "Full_Group_Imputed", ]

# ============================================================================
# SUMMARY STATISTICS (Perseus-style)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("PERSEUS-STYLE ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n")
cat("Imputation method:", imputation_method, "\n")
if (imputation_method == "perseus_downshift") {
  cat("  Downshift factor:", downshift_factor, "\n")
  cat("  Mean shift factor:", downshift_mean_shift, "\n")
}
cat("Total proteins analyzed:", nrow(results), "\n")
cat("Proteins with valid statistical tests:", nrow(complete_results), "\n")
cat("High quality results (≥", min_detection_threshold, " replicates both groups):", nrow(high_quality_results), "\n")
cat("Medium quality results (≥", min_detection_threshold, " replicates total):", nrow(medium_quality_results), "\n")
cat("Significant proteins (p <", significance_threshold, ", |FC| >", fold_change_threshold, "):", nrow(significant_results), "\n")
cat("  - Enriched in GNPT:", nrow(enriched_GNPT), "\n")
cat("  - Enriched in WT:", nrow(enriched_wt), "\n")

# Detailed imputation breakdown
cat("\nImputation Breakdown:\n")
cat("No imputation needed:", nrow(no_imputation), "\n")
cat("Partial imputation:", nrow(partial_imputation), "\n")
cat("Full group imputed:", nrow(full_group_imputed), "\n")

cat("\nData Quality Distribution:\n")
print(table(results$Data_Quality))

# ============================================================================
# PERSEUS-STYLE VISUALIZATION WITH TVAG LABELS
# ============================================================================

cat("\nCreating Perseus-style visualizations with TVAG labels...\n")

# Enhanced volcano plot with imputation categories and TVAG labels
if (nrow(complete_results) > 0) {
  # Add significance and enrichment categories
  complete_results$Category <- "Not Significant"
  
  # Upper right quadrant: positive fold change AND significant p-value = Enriched in GNPT
  complete_results$Category[
    !is.na(complete_results$P_Value) & 
      complete_results$P_Value < significance_threshold &
      complete_results$Log2_FC > fold_change_threshold
  ] <- "Enriched in GNPT"
  
  # Upper left quadrant: negative fold change AND significant p-value = Enriched in WT  
  complete_results$Category[
    !is.na(complete_results$P_Value) & 
      complete_results$P_Value < significance_threshold &
      complete_results$Log2_FC < -fold_change_threshold
  ] <- "Enriched in WT"
  
  # Enhanced categories combining enrichment and imputation
  complete_results$Enhanced_Category <- paste(
    complete_results$Category,
    complete_results$Imputation_Category,
    sep = "_"
  )
  
  # Create color palette for enhanced categories
  color_palette <- c(
    "Enriched in GNPT_No_Imputation" = "red",
    "Enriched in GNPT_Partial_Imputation" = "orange", 
    "Enriched in GNPT_Full_Group_Imputed" = "darkred",
    "Enriched in WT_No_Imputation" = "blue",
    "Enriched in WT_Partial_Imputation" = "lightblue",
    "Enriched in WT_Full_Group_Imputed" = "darkblue",
    "Not Significant_No_Imputation" = "lightgray",
    "Not Significant_Partial_Imputation" = "gray",
    "Not Significant_Full_Group_Imputed" = "darkgray"
  )
  
  # Create hover text
  complete_results$hover_text <- paste(
    "Gene:", complete_results$Protein_Gene,
    "<br>Description:", substr(complete_results$Protein_Description, 1, 60),
    "<br>Log2 FC:", round(complete_results$Log2_FC, 3),
    "<br>P-value:", format(complete_results$P_Value, scientific = TRUE, digits = 3),
    "<br>Adj P-value:", format(complete_results$Adj_P_Value, scientific = TRUE, digits = 3),
    "<br>Category:", complete_results$Category,
    "<br>Data Quality:", complete_results$Data_Quality,
    "<br>Imputation:", complete_results$Imputation_Category
  )
  
  # Identify proteins enriched in GNPT for labeling (avoid WT enrichment as requested)
  proteins_to_label <- complete_results[complete_results$Category == "Enriched in GNPT", ]
  
  # Static volcano plot for PNG/PDF export with TVAG labels
  static_volcano <- ggplot(complete_results, aes(x = Log2_FC, y = -log10(P_Value))) +
    geom_point(aes(color = Category), alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "Not Significant" = "gray", 
      "Enriched in GNPT" = "red", 
      "Enriched in WT" = "blue"
    )) +
    geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), 
               linetype = "dashed", color = "black", alpha = 0.7) +
    geom_hline(yintercept = -log10(significance_threshold), 
               linetype = "dashed", color = "black", alpha = 0.7) +
    # Add TVAG labels for proteins enriched in GNPT only
    geom_text_repel(
      data = proteins_to_label,
      aes(label = Protein_Gene),
      color = "darkred",
      size = 3,
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "darkred",
      segment.size = 0.3,
      segment.alpha = 0.6,
      min.segment.length = 0.1,
      seed = 123  # For reproducible label positioning
    ) +
    labs(
      title = "Perseus-style Volcano Plot: GNPT vs WT",
      subtitle = paste("Enriched GNPT proteins labeled with TVAG numbers | P-value <", significance_threshold, "| |Log2FC| >", fold_change_threshold),
      x = "Log2 Fold Change (GNPT/WT)",
      y = "-log10(P-value)",
      color = "Protein Category"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.margin = margin(20, 20, 20, 20)  # Add margins to accommodate labels
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
  # Add text annotations for quadrants (keeping the general quadrant labels)
  max_y <- max(-log10(complete_results$P_Value), na.rm = TRUE)
  max_x <- max(abs(complete_results$Log2_FC), na.rm = TRUE)
  
  static_volcano <- static_volcano +
    annotate("text", x = max_x * 0.7, y = max_y * 0.9, 
             label = "Enriched in\nGNPT", size = 4, color = "red", fontface = "bold") +
    annotate("text", x = -max_x * 0.7, y = max_y * 0.9, 
             label = "Enriched in\nWT", size = 4, color = "blue", fontface = "bold") +
    annotate("text", x = 0, y = max_y * 0.1, 
             label = "Not Significant", size = 4, color = "gray", fontface = "italic")
  
  print(static_volcano)
  
  # Create a version with enhanced labeling for very significant proteins
  # (Top proteins by combined significance and fold change)
  if (nrow(proteins_to_label) > 0) {
    # Calculate a combined score for ranking (p-value * fold change)
    proteins_to_label$combined_score <- -log10(proteins_to_label$P_Value) * abs(proteins_to_label$Log2_FC)
    
    # Alternative version: Label only top proteins to avoid overcrowding
    top_proteins <- proteins_to_label[order(proteins_to_label$combined_score, decreasing = TRUE), ]
    n_top_to_label <- min(20, nrow(top_proteins))  # Label top 20 or fewer if less available
    top_proteins_subset <- top_proteins[1:n_top_to_label, ]
    
    static_volcano_top <- ggplot(complete_results, aes(x = Log2_FC, y = -log10(P_Value))) +
      geom_point(aes(color = Category), alpha = 0.7, size = 2) +
      scale_color_manual(values = c(
        "Not Significant" = "gray", 
        "Enriched in GNPT" = "red", 
        "Enriched in WT" = "blue"
      )) +
      geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), 
                 linetype = "dashed", color = "black", alpha = 0.7) +
      geom_hline(yintercept = -log10(significance_threshold), 
                 linetype = "dashed", color = "black", alpha = 0.7) +
      # Add TVAG labels for top proteins enriched in GNPT only
      geom_text_repel(
        data = top_proteins_subset,
        aes(label = Protein_Gene),
        color = "darkred",
        size = 3,
        fontface = "bold",
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "darkred",
        segment.size = 0.3,
        segment.alpha = 0.6,
        min.segment.length = 0.1,
        seed = 123
      ) +
      labs(
        title = "Perseus-style Volcano Plot: GNPT vs WT (Top Enriched Labeled)",
        subtitle = paste("Top", n_top_to_label, "GNPT-enriched proteins labeled with TVAG numbers | P-value <", significance_threshold, "| |Log2FC| >", fold_change_threshold),
        x = "Log2 Fold Change (GNPT/WT)",
        y = "-log10(P-value)",
        color = "Protein Category"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(20, 20, 20, 20)
      ) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
      annotate("text", x = max_x * 0.7, y = max_y * 0.9, 
               label = "Enriched in\nGNPT", size = 4, color = "red", fontface = "bold") +
      annotate("text", x = -max_x * 0.7, y = max_y * 0.9, 
               label = "Enriched in\nWT", size = 4, color = "blue", fontface = "bold") +
      annotate("text", x = 0, y = max_y * 0.1, 
               label = "Not Significant", size = 4, color = "gray", fontface = "italic")
    
    print(static_volcano_top)
  }
  
  # Perseus-style interactive volcano plot (keeping original for web use)
  perseus_volcano <- plot_ly(
    data = complete_results,
    x = ~Log2_FC,
    y = ~-log10(P_Value),
    color = ~Enhanced_Category,
    colors = color_palette,
    text = ~hover_text,
    hovertemplate = "%{text}<extra></extra>",
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 6,
      opacity = 0.7,
      line = list(width = 0.5, color = "white")
    )
  ) %>%
    layout(
      title = list(
        text = paste0(
          "Perseus-style Volcano Plot: GNPT vs WT<br>",
          "<sub>Method: ", imputation_method, " | P-value threshold: ", significance_threshold, 
          " | FC threshold: ", fold_change_threshold, "</sub>"
        ),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Log2 Fold Change (GNPT/WT)",
        zeroline = TRUE,
        zerolinecolor = "black",
        zerolinewidth = 1
      ),
      yaxis = list(
        title = "-log10(P-value)",
        zeroline = FALSE
      ),
      hovermode = "closest",
      showlegend = TRUE,
      legend = list(
        orientation = "v",
        x = 1.02,
        y = 1
      ),
      plot_bgcolor = "white",
      shapes = list(
        # Threshold lines
        list(
          type = "line",
          x0 = fold_change_threshold, x1 = fold_change_threshold,
          y0 = 0, y1 = max(-log10(complete_results$P_Value), na.rm = TRUE),
          line = list(color = "black", width = 1, dash = "dash")
        ),
        list(
          type = "line",
          x0 = -fold_change_threshold, x1 = -fold_change_threshold,
          y0 = 0, y1 = max(-log10(complete_results$P_Value), na.rm = TRUE),
          line = list(color = "black", width = 1, dash = "dash")
        ),
        list(
          type = "line",
          x0 = min(complete_results$Log2_FC, na.rm = TRUE), 
          x1 = max(complete_results$Log2_FC, na.rm = TRUE),
          y0 = -log10(significance_threshold), 
          y1 = -log10(significance_threshold),
          line = list(color = "black", width = 1, dash = "dash")
        )
      )
    ) %>%
    config(
      toImageButtonOptions = list(
        format = "png",
        filename = "perseus_volcano_plot",
        height = 700,
        width = 1000,
        scale = 2
      )
    )
  
  print(perseus_volcano)
}

# Imputation quality assessment plot
imputation_plot <- ggplot(complete_results, aes(x = Imputation_Category, fill = Data_Quality)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  labs(
    title = "Imputation Usage by Data Quality",
    x = "Imputation Category",
    y = "Number of Proteins",
    fill = "Data Quality"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(imputation_plot)

# Distribution comparison before/after imputation
before_after_data <- data.frame(
  Values = c(as.vector(GNPT_filtered), as.vector(wt_filtered),
             as.vector(GNPT_imputed), as.vector(wt_imputed)),
  Type = rep(c("Before Imputation", "After Imputation"), each = length(c(as.vector(GNPT_filtered), as.vector(wt_filtered)))),
  Group = rep(rep(c("GNPT", "WT"), c(length(as.vector(GNPT_filtered)), length(as.vector(wt_filtered)))), 2)
)

# Remove NAs for plotting
before_after_data <- before_after_data[!is.na(before_after_data$Values), ]

distribution_plot <- ggplot(before_after_data, aes(x = Values, fill = Type)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~Group) +
  labs(
    title = "Distribution Before and After Perseus-style Imputation",
    x = "Log2 Intensity",
    y = "Density",
    fill = "Imputation Status"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Before Imputation" = "lightblue", "After Imputation" = "orange"))

print(distribution_plot)

# ============================================================================
# PERSEUS-STYLE EXPORT RESULTS WITH TVAG LABELS
# ============================================================================

cat("\nExporting Perseus-style results with TVAG-labeled plots...\n")

# Create output directory if it doesn't exist
if (!dir.exists("perseus_proteomics_results")) {
  dir.create("perseus_proteomics_results")
}

# Save all results with Perseus-style naming
write.csv(results, "perseus_proteomics_results/perseus_complete_analysis.csv", row.names = FALSE)
write.csv(complete_results, "perseus_proteomics_results/perseus_valid_statistical_results.csv", row.names = FALSE)
write.csv(high_quality_results, "perseus_proteomics_results/perseus_high_quality_results.csv", row.names = FALSE)
write.csv(medium_quality_results, "perseus_proteomics_results/perseus_medium_quality_results.csv", row.names = FALSE)
write.csv(significant_results, "perseus_proteomics_results/perseus_significant_proteins.csv", row.names = FALSE)
write.csv(enriched_GNPT, "perseus_proteomics_results/perseus_proteins_enriched_in_GNPT.csv", row.names = FALSE)
write.csv(enriched_wt, "perseus_proteomics_results/perseus_proteins_enriched_in_WT.csv", row.names = FALSE)

# Imputation-specific results
write.csv(no_imputation, "perseus_proteomics_results/perseus_no_imputation_proteins.csv", row.names = FALSE)
write.csv(partial_imputation, "perseus_proteomics_results/perseus_partial_imputation_proteins.csv", row.names = FALSE)
write.csv(full_group_imputed, "perseus_proteomics_results/perseus_full_group_imputed_proteins.csv", row.names = FALSE)

# Save plots with TVAG labels - PNG and PDF outputs
# All enriched GNPT proteins labeled
ggsave("perseus_proteomics_results/perseus_volcano_plot_TVAG_labeled.png", static_volcano, 
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave("perseus_proteomics_results/perseus_volcano_plot_TVAG_labeled.pdf", static_volcano, 
       width = 14, height = 10, device = "pdf")

# Top enriched proteins labeled (if available)
if (exists("static_volcano_top")) {
  ggsave("perseus_proteomics_results/perseus_volcano_plot_TOP_TVAG_labeled.png", static_volcano_top, 
         width = 14, height = 10, dpi = 300, bg = "white")
  ggsave("perseus_proteomics_results/perseus_volcano_plot_TOP_TVAG_labeled.pdf", static_volcano_top, 
         width = 14, height = 10, device = "pdf")
}

# Save other plots
ggsave("perseus_proteomics_results/perseus_imputation_quality_plot.png", imputation_plot, width = 10, height = 8, dpi = 300)
ggsave("perseus_proteomics_results/perseus_distribution_comparison.png", distribution_plot, width = 12, height = 8, dpi = 300)

# Save interactive plots
if (exists("perseus_volcano")) {
  saveWidget(perseus_volcano, "perseus_proteomics_results/perseus_interactive_volcano.html", selfcontained = TRUE)
}

# Create a summary of labeled proteins
if (exists("proteins_to_label") && nrow(proteins_to_label) > 0) {
  labeled_proteins_summary <- proteins_to_label %>%
    select(Protein_Gene, Protein_Description, Log2_FC, P_Value, Adj_P_Value, 
           Data_Quality, Imputation_Category) %>%
    arrange(desc(abs(Log2_FC))) %>%
    mutate(
      Log2_FC = round(Log2_FC, 3),
      P_Value = format(P_Value, scientific = TRUE, digits = 3),
      Adj_P_Value = format(Adj_P_Value, scientific = TRUE, digits = 3)
    )
  
  write.csv(labeled_proteins_summary, 
            "perseus_proteomics_results/TVAG_labeled_proteins_summary.csv", 
            row.names = FALSE)
  
  cat("Summary of TVAG-labeled proteins:\n")
  cat("Total proteins labeled:", nrow(proteins_to_label), "\n")
  cat("Proteins labeled in 'All' version:", nrow(proteins_to_label), "\n")
  if (exists("top_proteins_subset")) {
    cat("Proteins labeled in 'Top' version:", nrow(top_proteins_subset), "\n")
  }
}

cat("\nTVAG-labeled volcano plots saved:\n")
cat("- perseus_volcano_plot_TVAG_labeled.png (all enriched GNPT proteins labeled)\n")
cat("- perseus_volcano_plot_TVAG_labeled.pdf (all enriched GNPT proteins labeled)\n")
if (exists("static_volcano_top")) {
  cat("- perseus_volcano_plot_TOP_TVAG_labeled.png (top enriched proteins labeled)\n")
  cat("- perseus_volcano_plot_TOP_TVAG_labeled.pdf (top enriched proteins labeled)\n")
}
cat("- TVAG_labeled_proteins_summary.csv (summary of labeled proteins)\n")

# ============================================================================
# PERSEUS-STYLE DETAILED SUMMARY REPORT
# ============================================================================

# Create comprehensive Perseus-style summary report
perseus_summary <- paste(
  "PERSEUS-STYLE PROTEOMICS ANALYSIS SUMMARY REPORT",
  paste(rep("=", 70), collapse = ""),
  "",
  paste("Analysis Date:", Sys.Date()),
  paste("Input File:", file_path),
  paste("Imputation Method:", imputation_method),
  "",
  "IMPUTATION PARAMETERS:",
  if (imputation_method == "perseus_downshift") {
    paste(
      paste("  Width Factor (dispersion):", downshift_factor),
      paste("  Mean Shift Factor:", downshift_mean_shift),
      paste("  Strategy: Downshifted normal distribution"),
      paste("  Description: Missing values replaced by random sampling from"),
      paste("    a normal distribution with reduced mean and variance"),
      sep = "\n"
    )
  } else if (imputation_method == "perseus_normal") {
    paste(
      paste("  Strategy: Normal distribution from valid values"),
      paste("  Description: Missing values replaced by random sampling from"),
      paste("    the distribution of valid intensity values"),
      sep = "\n"
    )
  } else {
    paste("  Strategy:", imputation_method)
  },
  "",
  "ANALYSIS THRESHOLDS:",
  paste("  Significance threshold (p-value):", significance_threshold),
  paste("  Fold change threshold (log2):", fold_change_threshold),
  paste("  Minimum detection threshold:", min_detection_threshold),
  paste("  Valid values threshold:", valid_values_threshold),
  "",
  "DATA SUMMARY:",
  paste("  Total proteins in dataset:", nrow(data)),
  paste("  Proteins with sufficient data:", nrow(results)),
  paste("  Proteins excluded (insufficient data):", nrow(data) - nrow(results)),
  "",
  "STATISTICAL RESULTS:",
  paste("  Proteins with valid statistical tests:", nrow(complete_results)),
  paste("  High quality results (≥", min_detection_threshold, " replicates both groups):", nrow(high_quality_results)),
  paste("  Medium quality results (≥", min_detection_threshold, " replicates total):", nrow(medium_quality_results)),
  paste("  Low quality results:", nrow(complete_results) - nrow(high_quality_results) - nrow(medium_quality_results)),
  "",
  "SIGNIFICANCE ANALYSIS:",
  paste("  Significant proteins (raw p-value <", significance_threshold, "):", nrow(significant_results)),
  paste("    - Enriched in GNPT (upper right quadrant):", nrow(enriched_GNPT)),
  paste("    - Enriched in WT (upper left quadrant):", nrow(enriched_wt)),
  paste("  Non-significant proteins:", nrow(complete_results) - nrow(significant_results)),
  "",
  "IMPUTATION BREAKDOWN:",
  paste("  No imputation needed:", nrow(no_imputation)),
  paste("  Partial imputation:", nrow(partial_imputation)),
  paste("  Full group imputed:", nrow(full_group_imputed)),
  paste("  Total imputed values:", sum(GNPT_imputation_mask) + sum(wt_imputation_mask)),
  "",
  "PERSEUS-STYLE DATA QUALITY CATEGORIES:",
  paste("  High Quality: Both groups have ≥", min_detection_threshold, "valid replicates"),
  paste("  Medium Quality: Combined groups have ≥", min_detection_threshold, "valid replicates"),
  paste("  Low Quality: Insufficient valid replicates, relies heavily on imputation"),
  "",
  "IMPUTATION PHILOSOPHY (Perseus-style):",
  "  - Missing values are assumed to be missing not at random (MNAR)",
  "  - Imputation uses downshifted distributions to simulate low-abundance proteins",
  "  - Statistical tests preferentially use original data when available",
  "  - Quality metrics help interpret reliability of results",
  "",
  "FILES GENERATED:",
  "  Analysis Results:",
  "    - perseus_complete_analysis.csv (all proteins)",
  "    - perseus_valid_statistical_results.csv (proteins with p-values)",
  "    - perseus_high_quality_results.csv (high confidence results)",
  "    - perseus_medium_quality_results.csv (medium confidence results)",
  "    - perseus_significant_proteins.csv (significant results)",
  "    - perseus_proteins_enriched_in_GNPT.csv",
  "    - perseus_proteins_enriched_in_WT.csv",
  "",
  "  Imputation Analysis:",
  "    - perseus_no_imputation_proteins.csv",
  "    - perseus_partial_imputation_proteins.csv", 
  "    - perseus_full_group_imputed_proteins.csv",
  "",
  "  Visualizations with TVAG Labels:",
  "    - perseus_volcano_plot_TVAG_labeled.png (all enriched GNPT proteins labeled)",
  "    - perseus_volcano_plot_TVAG_labeled.pdf (publication-ready with all labels)",
  "    - perseus_volcano_plot_TOP_TVAG_labeled.png (top enriched proteins only)",
  "    - perseus_volcano_plot_TOP_TVAG_labeled.pdf (cleaner version with top proteins)",
  "    - perseus_interactive_volcano.html (interactive volcano plot)",
  "    - perseus_imputation_quality_plot.png",
  "    - perseus_distribution_comparison.png",
  "    - TVAG_labeled_proteins_summary.csv (summary of all labeled proteins)",
  "",
  "VOLCANO PLOT INTERPRETATION:",
  "  - Upper Right Quadrant: Proteins enriched in GNPT (positive fold change + significant p-value)",
  "  - Upper Left Quadrant: Proteins enriched in WT (negative fold change + significant p-value)",  
  "  - Lower Quadrants: Non-significant proteins regardless of fold change",
  "  - Threshold lines: Vertical = fold change cutoff, Horizontal = p-value cutoff",
  "  - TVAG Labels: Only proteins enriched in GNPT are labeled with TVAG numbers",
  "",
  "TVAG LABELING STRATEGY:",
  "  - All Version: Labels all significantly enriched GNPT proteins",
  "  - Top Version: Labels top 20 most significant enriched GNPT proteins",
  "  - Ranking: Based on combined score (-log10(p-value) × |fold change|)",
  "  - No WT labeling: As requested, WT-enriched proteins are not labeled",
  "",
  "QUALITY INTERPRETATION GUIDELINES:",
  "",
  "1. HIGH QUALITY RESULTS (recommended for publication):",
  "   - Both groups have sufficient original data",
  "   - Statistical tests based on measured values",
  "   - Most reliable for biological interpretation",
  "",
  "2. MEDIUM QUALITY RESULTS (interpret with caution):",
  "   - Some imputation used but sufficient total data",
  "   - May be biologically meaningful but require validation",
  "   - Consider fold change magnitude and biological context",
  "",
  "3. LOW QUALITY RESULTS (exploratory only):",
  "   - Heavy reliance on imputation",
  "   - Useful for hypothesis generation",
  "   - Require experimental validation",
  "",
  "4. IMPUTATION CATEGORIES:",
  "   - No_Imputation: Most reliable, no missing values",
  "   - Partial_Imputation: Some missing values, mixed reliability", 
  "   - Full_Group_Imputed: One group entirely imputed, lowest confidence",
  "",
  "RECOMMENDED NEXT STEPS:",
  "1. Focus on high-quality significant results for detailed analysis",
  "2. Validate medium-quality results with independent experiments",
  "3. Consider pathway analysis with significant protein lists",
  "4. Examine imputation distribution plots for quality assessment",
  "5. Use TVAG-labeled plots for presentation and publication",
  "",
  paste("Analysis completed:", Sys.time()),
  "",
  sep = "\n"
)

writeLines(perseus_summary, "perseus_proteomics_results/perseus_analysis_summary.txt")

# ============================================================================
# ADDITIONAL PERSEUS-STYLE QUALITY METRICS
# ============================================================================

cat("\nCalculating additional Perseus-style quality metrics...\n")

# Calculate imputation statistics
imputation_stats <- data.frame(
  Metric = c(
    "Total proteins analyzed",
    "Proteins with complete data (no imputation)",
    "Proteins with partial imputation",
    "Proteins with full group imputation",
    "Total values imputed",
    "Percentage of values imputed",
    "GNPT values imputed",
    "WT values imputed"
  ),
  Value = c(
    nrow(results),
    nrow(no_imputation),
    nrow(partial_imputation),
    nrow(full_group_imputed),
    sum(GNPT_imputation_mask) + sum(wt_imputation_mask),
    round(100 * (sum(GNPT_imputation_mask) + sum(wt_imputation_mask)) / 
            (length(GNPT_imputation_mask) + length(wt_imputation_mask)), 2),
    sum(GNPT_imputation_mask),
    sum(wt_imputation_mask)
  )
)

write.csv(imputation_stats, "perseus_proteomics_results/perseus_imputation_statistics.csv", row.names = FALSE)

# Quality distribution by fold change
quality_fc_breakdown <- complete_results %>%
  mutate(
    FC_Category = case_when(
      abs(Log2_FC) < 0.5 ~ "Low FC (<0.5)",
      abs(Log2_FC) < 1.0 ~ "Medium FC (0.5-1.0)",
      abs(Log2_FC) < 2.0 ~ "High FC (1.0-2.0)",
      TRUE ~ "Very High FC (>2.0)"
    ),
    Significance = ifelse(P_Value < significance_threshold, "Significant", "Not Significant")
  ) %>%
  count(Data_Quality, FC_Category, Significance) %>%
  pivot_wider(names_from = c(FC_Category, Significance), values_from = n, values_fill = 0)

write.csv(quality_fc_breakdown, "perseus_proteomics_results/perseus_quality_fc_breakdown.csv", row.names = FALSE)

cat("\nPereus-style analysis with TVAG labeling completed successfully!\n")
cat("Results saved in 'perseus_proteomics_results' folder\n")
cat("Check 'perseus_analysis_summary.txt' for detailed interpretation guidelines\n")
cat("\nKey Perseus-style features implemented:\n")
cat("✓ Downshifted normal distribution imputation\n")
cat("✓ Quality-based statistical analysis\n")
cat("✓ Imputation category tracking\n")
cat("✓ Raw p-value volcano plots\n")
cat("✓ TVAG labels for enriched GNPT proteins\n")
cat("✓ Two versions: All labeled vs Top 20 labeled\n")
cat("✓ Comprehensive quality metrics\n")
cat("✓ Perseus-style interpretation guidelines\n")
cat("\nNew TVAG labeling features:\n")
cat("✓ Labels only GNPT-enriched proteins (no WT labeling)\n")
cat("✓ Anti-overlapping labels using ggrepel\n")
cat("✓ Two plot versions for different needs\n")
cat("✓ Summary file of all labeled proteins\n")