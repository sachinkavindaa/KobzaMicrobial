# Load required libraries
library(phyloseq)
library(DESeq2)
library(pheatmap)

# Load your cleaned phyloseq object
ps_mock_analyze <- readRDS("~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")

# Ensure TRT factor is properly ordered
sample_data(ps_mock_analyze)$TRT <- factor(sample_data(ps_mock_analyze)$TRT,
                                           levels = c("Control", "COL", "TL", "Tallow", "Corn Oil"))

# Convert phyloseq object to DESeq2 object
dds <- phyloseq_to_deseq2(ps_mock_analyze, ~ TRT)

# Estimate size factors using geometric means
geoMeans <- apply(counts(dds), 1, function(x) exp(sum(log(x[x > 0])) / length(x)))
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)

# Run DESeq2
dds <- DESeq(dds, fitType = "local")

# Apply Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Apply row-wise Z-score transformation
vsd_mat_z <- t(scale(t(vsd_mat)))

# Identify significant ASVs (adjusted p-value < 0.05)
res_all <- results(dds)
alpha <- 0.05
sig_asv_names <- rownames(res_all)[which(res_all$padj < alpha)]

# Subset matrix to only significant ASVs
vsd_mat_sig <- vsd_mat_z[sig_asv_names, ]

# Taxonomic labeling (use Family if Genus is missing)
taxa <- as.data.frame(tax_table(ps_mock_analyze))
taxa$Genus <- as.character(taxa$Genus)
taxa$Family <- as.character(taxa$Family)

# Replace NA or blank Genus with Family
missing_genus <- is.na(taxa$Genus) | taxa$Genus == ""
taxa$Genus[missing_genus] <- taxa$Family[missing_genus]

# Replace remaining NAs with "Unclassified"
taxa$Genus[is.na(taxa$Genus) | taxa$Genus == ""] <- "Unclassified"

# Create informative row labels: Genus_ASVxxxx
genus_labels <- paste0(taxa[sig_asv_names, "Genus"], "_", sig_asv_names)
rownames(vsd_mat_sig) <- genus_labels

# Order samples by treatment
sample_order <- order(colData(vsd)[, "TRT"])
vsd_mat_sig <- vsd_mat_sig[, sample_order]
annotation_col <- as.data.frame(colData(vsd)[sample_order, "TRT", drop = FALSE])

# Plot heatmap
pheatmap(vsd_mat_sig,
         annotation_col = annotation_col,
         fontsize = 8,
         fontsize_row = 5.5,
         show_rownames = TRUE,
         cluster_cols = FALSE,  # keep sample order by TRT
         cluster_rows = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("Low", "Average", "High"),
         main = "Differentially Abundant ASVs Across Treatments"
)


##

