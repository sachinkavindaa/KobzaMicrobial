# Load necessary libraries
library(phyloseq)
library(DESeq2)
library(UpSetR)
library(ggVennDiagram)


# Load your phyloseq object
ps_mock_analyze <- readRDS("~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")

# Define treatment groups
treatments <- c("Control", "COL", "TL", "Tallow", "Corn Oil")

# Initialize a list to hold significant ASVs from each comparison
sig_asv_list <- list()

# Loop through all pairwise treatment combinations
for (i in 1:(length(treatments) - 1)) {
  for (j in (i + 1):length(treatments)) {
    
    # Extract treatment names
    trt1 <- treatments[i]
    trt2 <- treatments[j]
    comp_name <- paste(trt2, "vs", trt1, sep = "_")
    
    cat("Running:", comp_name, "\n")
    
    # Subset phyloseq object
    ps_sub <- subset_samples(ps_mock_analyze, TRT %in% c(trt1, trt2))
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

    # Create DESeq2 object
    dds <- phyloseq_to_deseq2(ps_sub, ~TRT)
    
    # Estimate size factors
    geoMeans <- apply(counts(dds), 1, function(x) exp(sum(log(x[x > 0])) / length(x)))
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    
    # Run DESeq2
    dds <- DESeq(dds, fitType = "local")
    
    # Extract results with contrast
    res <- results(dds, contrast = c("TRT", trt2, trt1), cooksCutoff = FALSE)
    
    # Filter significant ASVs
    sig_asv <- rownames(res)[which(res$padj < 0.05)]
    
    # Store in list
    sig_asv_list[[comp_name]] <- sig_asv
  }
}

# Convert the list to a binary matrix for UpSet plot
asv_names <- unique(unlist(sig_asv_list))
binary_matrix <- sapply(sig_asv_list, function(x) as.integer(asv_names %in% x))
rownames(binary_matrix) <- asv_names

# Convert to data.frame
upset_data <- as.data.frame(binary_matrix)

# Plot UpSet
upset(upset_data, 
      nsets = length(sig_asv_list), 
      nintersects = 20,
      order.by = "freq", 
      main.bar.color = "steelblue", 
      sets.bar.color = "#3b5998",
      point.size = 2.0,
      line.size = 1.0,
      text.scale = c(1.4, 1.5, 1.4, 1.2, 1.0, 1.4),
      sets = rev(names(sig_asv_list)))

selected_comparisons <- c(
  "Corn Oil_vs_Control",
  "Tallow_vs_Control",
  "Corn Oil_vs_TL",
  "Corn Oil_vs_COL",
  "Tallow_vs_TL"
)

venn_list <- sig_asv_list[selected_comparisons]

# Convert binary_matrix to a data.frame and add ASV IDs
upset_df <- as.data.frame(binary_matrix)
upset_df$ASV <- rownames(binary_matrix)

# Function to extract ASVs for each unique intersection
get_upset_intersections <- function(data) {
  sets <- colnames(data)[colnames(data) != "ASV"]
  
  # Create all unique intersection patterns (columns with 1's)
  intersections <- data %>%
    pivot_longer(-ASV, names_to = "Comparison", values_to = "Present") %>%
    group_by(ASV) %>%
    summarise(Pattern = paste(Comparison[Present == 1], collapse = " & ")) %>%
    group_by(Pattern) %>%
    summarise(ASVs = list(ASV), Count = n()) %>%
    arrange(desc(Count))
  
  return(intersections)
}

# Run function
intersection_results <- get_upset_intersections(upset_df)

# View summary table
intersection_results


tax_tab <- as.data.frame(tax_table(ps_mock_analyze))

# Add ASV IDs as a column (rownames are sequence IDs)
tax_tab$ASV <- rownames(tax_tab)

# Inspect a few rows
head(tax_tab)

upset_df <- as.data.frame(binary_matrix)
upset_df$ASV <- rownames(binary_matrix)

get_upset_intersections <- function(data) {
  data %>%
    pivot_longer(-ASV, names_to = "Comparison", values_to = "Present") %>%
    group_by(ASV) %>%
    summarise(Pattern = if_else(sum(Present) > 0,
                                paste(sort(Comparison[Present == 1]), collapse = " & "), NA_character_)) %>%
    drop_na() %>%
    group_by(Pattern) %>%
    summarise(ASVs = list(ASV), Count = n())
}

intersection_results <- get_upset_intersections(upset_df)

intersection_tax <- intersection_results %>%
  unnest_longer(ASVs) %>%
  rename(ASV = ASVs) %>%
  left_join(tax_tab, by = "ASV")

head(intersection_tax)

write.csv(intersection_tax,
          file = "~/Documents/Kobza/intersection_taxonomy_per_intersection.csv",
          row.names = FALSE)


