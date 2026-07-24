# Radial Co-IP Interactome Network Visualization
# Creates the specific radial layout shown in the image
# Author: Generated for GNPT, MPR3, MPR7 bait proteins
# Date: 2025

# Load required libraries
library(igraph)
library(ggraph)
library(tidyverse)
library(RColorBrewer)

# Try to load Cairo for high-quality output
cairo_available <- FALSE
tryCatch({
  library(Cairo)
  cairo_available <- TRUE
  options(bitmapType='cairo')
  cat("Cairo loaded successfully for high-quality graphics\n")
}, error = function(e) {
  cat("Cairo not available, using standard graphics\n")
})

# Function to read and process co-IP data (same as before)
read_coip_data <- function(file_path, bait_name, bait_id) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat("Successfully read", file_path, "- Rows:", nrow(data), "\n")
  
  # Filter for significant interactions
  significant <- data %>%
    filter(!is.na(Log2_FC), Log2_FC > 0.58) %>%
    arrange(desc(Log2_FC))
  
  cat("Filtered to", nrow(significant), "significant interactions for", bait_name, "\n")
  
  significant$Bait <- bait_name
  significant$Bait_ID <- bait_id
  return(significant)
}

# Check for data files
data_files <- c("GNPT_coIP.csv", "MPR3_coIP.csv", "MPR7_coIP.csv")
missing_files <- data_files[!file.exists(data_files)]

if (length(missing_files) > 0) {
  stop(paste("Missing data files:", paste(missing_files, collapse = ", ")))
}

# Read data
gnpt_data <- read_coip_data("GNPT_coIP.csv", "GNPT", "TVAG_166090")
mpr3_data <- read_coip_data("MPR3_coIP.csv", "MPR3", "TVAG_498650")
mpr7_data <- read_coip_data("MPR7_coIP.csv", "MPR7", "TVAG_255730")

# Combine all data
all_data <- bind_rows(gnpt_data, mpr3_data, mpr7_data)

# Create edge list
create_edge_list <- function(data) {
  protein_col <- "Protein_Gene"
  if (!"Protein_Gene" %in% names(data)) {
    protein_col <- names(data)[grep("Gene|Protein", names(data), ignore.case = TRUE)[1]]
    if (is.na(protein_col)) {
      for (col in names(data)) {
        if (any(grepl("TVAG_", data[[col]], na.rm = TRUE))) {
          protein_col <- col
          break
        }
      }
    }
  }
  
  edges <- data.frame(
    from = data$Bait_ID,
    to = data[[protein_col]],
    weight = data$Log2_FC,
    bait = data$Bait,
    stringsAsFactors = FALSE
  )
  
  edges <- edges[!is.na(edges$to) & edges$to != "", ]
  return(edges)
}

edges <- create_edge_list(all_data)

# Find overlapping proteins
gnpt_proteins <- unique(edges$to[edges$bait == "GNPT"])
mpr3_proteins <- unique(edges$to[edges$bait == "MPR3"])
mpr7_proteins <- unique(edges$to[edges$bait == "MPR7"])

three_way_overlap <- intersect(intersect(gnpt_proteins, mpr3_proteins), mpr7_proteins)
gnpt_mpr3_only <- setdiff(intersect(gnpt_proteins, mpr3_proteins), mpr7_proteins)
gnpt_mpr7_only <- setdiff(intersect(gnpt_proteins, mpr7_proteins), mpr3_proteins)
mpr3_mpr7_only <- setdiff(intersect(mpr3_proteins, mpr7_proteins), gnpt_proteins)

# Unique proteins for each bait
gnpt_unique <- setdiff(setdiff(gnpt_proteins, mpr3_proteins), mpr7_proteins)
mpr3_unique <- setdiff(setdiff(mpr3_proteins, gnpt_proteins), mpr7_proteins)
mpr7_unique <- setdiff(setdiff(mpr7_proteins, gnpt_proteins), mpr3_proteins)

cat("Overlap analysis:\n")
cat("Three-way overlap:", length(three_way_overlap), "\n")
cat("GNPT unique:", length(gnpt_unique), "\n")
cat("MPR3 unique:", length(mpr3_unique), "\n") 
cat("MPR7 unique:", length(mpr7_unique), "\n")

# Add overlap type to edges
edges$overlap_type <- case_when(
  edges$to %in% three_way_overlap ~ "three_way",
  edges$to %in% gnpt_mpr3_only ~ "gnpt_mpr3",
  edges$to %in% gnpt_mpr7_only ~ "gnpt_mpr7",
  edges$to %in% mpr3_mpr7_only ~ "mpr3_mpr7",
  TRUE ~ "unique"
)

# Create nodes
nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  stringsAsFactors = FALSE
)

nodes$type <- ifelse(nodes$name %in% c("TVAG_166090", "TVAG_498650", "TVAG_255730"), "bait", "prey")
nodes$bait_protein <- case_when(
  nodes$name == "TVAG_166090" ~ "GNPT",
  nodes$name == "TVAG_498650" ~ "MPR3",
  nodes$name == "TVAG_255730" ~ "MPR7",
  TRUE ~ "prey"
)

# Create graph
g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

# Set node colors matching the image
V(g)$color <- case_when(
  V(g)$name == "TVAG_166090" ~ "#E41A1C",      # GNPT - red
  V(g)$name == "TVAG_498650" ~ "#4DAF4A",      # MPR3 - green  
  V(g)$name == "TVAG_255730" ~ "#984EA3",      # MPR7 - purple
  V(g)$name %in% three_way_overlap ~ "#0066CC", # Three-way - blue
  V(g)$name %in% c(gnpt_mpr3_only, gnpt_mpr7_only, mpr3_mpr7_only) ~ "#87CEEB", # Two-way - light blue
  TRUE ~ "#D3D3D3"                             # Unique - gray
)

# Set edge colors based on bait source
E(g)$color <- case_when(
  E(g)$bait == "GNPT" ~ "#FF6B6B",     # Light red for GNPT edges
  E(g)$bait == "MPR3" ~ "#90EE90",     # Light green for MPR3 edges
  E(g)$bait == "MPR7" ~ "#DDA0DD",     # Light purple for MPR7 edges
  TRUE ~ "#CCCCCC"
)

# Set node sizes
V(g)$size <- case_when(
  V(g)$type == "bait" ~ 12,
  V(g)$name %in% three_way_overlap ~ 8,
  V(g)$name %in% c(gnpt_mpr3_only, gnpt_mpr7_only, mpr3_mpr7_only) ~ 6,
  TRUE ~ 4
)

# CUSTOM RADIAL LAYOUT FUNCTION
create_radial_layout <- function(graph) {
  n_nodes <- vcount(graph)
  layout_matrix <- matrix(0, nrow = n_nodes, ncol = 2)
  
  # Define bait positions in a triangle formation
  bait_positions <- list(
    "TVAG_166090" = c(-2, 2),    # GNPT - top left
    "TVAG_498650" = c(0, -2.5),  # MPR3 - bottom 
    "TVAG_255730" = c(2, 2)      # MPR7 - top right
  )
  
  # Position three-way overlap proteins in the center
  center_proteins <- V(graph)$name[V(graph)$name %in% three_way_overlap]
  if(length(center_proteins) > 0) {
    # Create a tight cluster in the center
    n_center <- length(center_proteins)
    center_angles <- seq(0, 2*pi, length.out = n_center + 1)[1:n_center]
    center_radius <- 0.5
    
    for(i in 1:length(center_proteins)) {
      node_idx <- which(V(graph)$name == center_proteins[i])
      layout_matrix[node_idx, 1] <- center_radius * cos(center_angles[i])
      layout_matrix[node_idx, 2] <- center_radius * sin(center_angles[i])
    }
  }
  
  # Position bait proteins
  for(bait_id in names(bait_positions)) {
    node_idx <- which(V(graph)$name == bait_id)
    if(length(node_idx) > 0) {
      layout_matrix[node_idx, ] <- bait_positions[[bait_id]]
    }
  }
  
  # Position unique proteins for each bait in radial pattern
  bait_data <- list(
    list(bait_id = "TVAG_166090", unique_proteins = gnpt_unique, center = bait_positions[["TVAG_166090"]], 
         angle_start = pi, angle_end = 0),  # Left semicircle
    list(bait_id = "TVAG_498650", unique_proteins = mpr3_unique, center = bait_positions[["TVAG_498650"]], 
         angle_start = pi/3, angle_end = 2*pi/3),  # Bottom area
    list(bait_id = "TVAG_255730", unique_proteins = mpr7_unique, center = bait_positions[["TVAG_255730"]], 
         angle_start = 0, angle_end = pi)  # Right semicircle
  )
  
  for(bait_info in bait_data) {
    unique_prots <- bait_info$unique_proteins
    if(length(unique_prots) > 0) {
      n_unique <- length(unique_prots)
      
      # Create multiple rings for better distribution
      rings <- ceiling(n_unique / 30)  # ~30 proteins per ring
      for(ring in 1:rings) {
        ring_proteins <- unique_prots[((ring-1)*30 + 1):min(ring*30, n_unique)]
        ring_proteins <- ring_proteins[!is.na(ring_proteins)]
        
        if(length(ring_proteins) > 0) {
          ring_radius <- 3 + ring * 1.5  # Increasing radius for each ring
          angles <- seq(bait_info$angle_start, bait_info$angle_end, 
                        length.out = length(ring_proteins) + 1)[1:length(ring_proteins)]
          
          for(i in 1:length(ring_proteins)) {
            node_idx <- which(V(graph)$name == ring_proteins[i])
            if(length(node_idx) > 0) {
              layout_matrix[node_idx, 1] <- bait_info$center[1] + ring_radius * cos(angles[i])
              layout_matrix[node_idx, 2] <- bait_info$center[2] + ring_radius * sin(angles[i])
            }
          }
        }
      }
    }
  }
  
  # Position two-way overlap proteins between their respective baits
  two_way_groups <- list(
    list(proteins = gnpt_mpr3_only, bait1 = "TVAG_166090", bait2 = "TVAG_498650"),
    list(proteins = gnpt_mpr7_only, bait1 = "TVAG_166090", bait2 = "TVAG_255730"), 
    list(proteins = mpr3_mpr7_only, bait1 = "TVAG_498650", bait2 = "TVAG_255730")
  )
  
  for(group in two_way_groups) {
    if(length(group$proteins) > 0) {
      pos1 <- bait_positions[[group$bait1]]
      pos2 <- bait_positions[[group$bait2]]
      mid_point <- c((pos1[1] + pos2[1])/2, (pos1[2] + pos2[2])/2)
      
      # Spread proteins around the midpoint
      n_proteins <- length(group$proteins)
      angles <- seq(0, 2*pi, length.out = n_proteins + 1)[1:n_proteins]
      radius <- 1.2
      
      for(i in 1:n_proteins) {
        node_idx <- which(V(graph)$name == group$proteins[i])
        if(length(node_idx) > 0) {
          layout_matrix[node_idx, 1] <- mid_point[1] + radius * cos(angles[i])
          layout_matrix[node_idx, 2] <- mid_point[2] + radius * sin(angles[i])
        }
      }
    }
  }
  
  return(layout_matrix)
}

# Create the radial layout
radial_layout <- create_radial_layout(g)

# Function to create high-quality plots
create_radial_plot <- function(filename, width = 16, height = 16, res = 600, show_labels = FALSE) {
  
  # PDF version
  if(grepl("\\.pdf$", filename)) {
    if(cairo_available) {
      cairo_pdf(filename, width = width, height = height, pointsize = 10)
    } else {
      pdf(filename, width = width, height = height, pointsize = 10)
    }
  } else {
    # PNG version
    if(cairo_available) {
      png(filename, width = width, height = height, units = "in", res = res, type = "cairo-png")
    } else {
      png(filename, width = width, height = height, units = "in", res = res)
    }
  }
  
  # Set up plotting parameters
  par(mar = c(1, 1, 3, 1), bg = "white")
  
  # Create the plot
  plot(g, 
       layout = radial_layout,
       vertex.size = V(g)$size,
       vertex.color = V(g)$color,
       vertex.frame.color = "black",
       vertex.frame.width = 0.5,
       vertex.label = if(show_labels) {
         case_when(
           V(g)$name == "TVAG_166090" ~ "GNPT",
           V(g)$name == "TVAG_498650" ~ "MPR3", 
           V(g)$name == "TVAG_255730" ~ "MPR7",
           TRUE ~ ""
         )
       } else {
         ""
       },
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       vertex.label.font = 2,
       vertex.label.dist = 0.5,
       edge.color = E(g)$color,
       edge.width = 0.5,
       edge.alpha = 0.7,
       edge.arrow.size = 0.3,
       edge.arrow.width = 0.8,
       edge.curved = 0.1,
       main = "Co-IP Interactome Network - Radial Layout",
       cex.main = 2.0,
       col.main = "#2E2E2E")
  
  dev.off()
  cat("Created:", filename, "\n")
}

# Create multiple versions of the radial plot
cat("Creating radial network visualizations...\n")

# High-resolution PDF version (main output)
create_radial_plot("coip_radial_network_HD.pdf", width = 20, height = 20, show_labels = FALSE)

# High-resolution PNG version
create_radial_plot("coip_radial_network_HD.png", width = 20, height = 20, res = 600, show_labels = FALSE)

# Labeled version for reference
create_radial_plot("coip_radial_network_labeled.pdf", width = 20, height = 20, show_labels = TRUE)

# Extra high resolution version
create_radial_plot("coip_radial_network_ultra_HD.pdf", width = 24, height = 24, show_labels = FALSE)

# GGRAPH VERSION for even better quality
tryCatch({
  library(ggplot2)
  library(ggraph)
  
  # Convert layout to data frame for ggraph
  layout_df <- data.frame(
    name = V(g)$name,
    x = radial_layout[, 1],
    y = radial_layout[, 2]
  )
  
  # Create tidy graph with layout
  g_tidy <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    left_join(layout_df, by = "name")
  
  # Create ggraph plot
  p_radial <- ggraph(g_tidy, x = x, y = y) +
    geom_edge_link(aes(color = bait), alpha = 0.6, width = 0.3,
                   arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                   end_cap = circle(2, 'mm')) +
    geom_node_point(aes(color = bait_protein, size = type), stroke = 0.3, alpha = 0.9) +
    geom_node_text(aes(label = case_when(
      name == "TVAG_166090" ~ "GNPT",
      name == "TVAG_498650" ~ "MPR3", 
      name == "TVAG_255730" ~ "MPR7",
      TRUE ~ ""
    )), size = 4, fontface = "bold", color = "black", nudge_y = 0.15) +
    scale_color_manual(values = c("GNPT" = "#E41A1C", "MPR3" = "#4DAF4A", "MPR7" = "#984EA3", "prey" = "#D3D3D3")) +
    scale_edge_color_manual(values = c("GNPT" = "#FF6B6B", "MPR3" = "#90EE90", "MPR7" = "#DDA0DD")) +
    scale_size_manual(values = c("bait" = 6, "prey" = 2)) +
    theme_graph(background = "white", base_family = "sans") +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5, color = "#2E2E2E")) +
    labs(title = "Co-IP Interactome Network - Radial Layout (ggraph)") +
    coord_fixed()
  
  # Save ggraph version
  ggsave("coip_radial_network_ggraph.pdf", p_radial, width = 16, height = 16, device = "pdf")
  ggsave("coip_radial_network_ggraph.png", p_radial, width = 16, height = 16, dpi = 600, device = "png")
  
  cat("Created ggraph versions\n")
  
}, error = function(e) {
  cat("Error creating ggraph version:", e$message, "\n")
})

# Summary
cat("\n=== Radial Network Visualization Complete ===\n")
cat("Files created:\n")
files_to_check <- c(
  "coip_radial_network_HD.pdf",
  "coip_radial_network_HD.png", 
  "coip_radial_network_labeled.pdf",
  "coip_radial_network_ultra_HD.pdf",
  "coip_radial_network_ggraph.pdf",
  "coip_radial_network_ggraph.png"
)

for(file in files_to_check) {
  if(file.exists(file)) {
    size_mb <- round(file.info(file)$size / 1024 / 1024, 2)
    cat("✓", file, paste0("(", size_mb, " MB)"), "\n")
  } else {
    cat("✗", file, "- NOT CREATED\n")
  }
}

cat("\nNetwork Statistics:\n")
cat("- Three bait proteins positioned as central hubs\n")
cat("- Three-way overlap proteins (", length(three_way_overlap), ") clustered in center\n")
cat("- Unique proteins radiate outward from their respective baits\n") 
cat("- Two-way overlaps positioned between corresponding baits\n")
cat("- Total nodes:", vcount(g), "| Total edges:", ecount(g), "\n")

cat("\n✨ RECOMMENDED OUTPUT: coip_radial_network_HD.pdf\n")
cat("This matches your requested radial layout with high resolution suitable for publication.\n")