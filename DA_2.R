# Load libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)

# Load your phyloseq object
load("~/Documents/Kobza/250/ps_mock_analyze_new.rds")
ps_mock_analyze
saveRDS(ps_mock_analyze, file = "~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")
ps_mock_analyze <- readRDS("~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")



# Ensure factor levels are correct
sample_data(ps_mock_analyze)$TRT <- factor(sample_data(ps_mock_analyze)$TRT, 
                                           levels = c("Control", "TL"))

# Remove NA values
ps_mock_analyze <- subset_samples(ps_mock_analyze, !is.na(TRT))

# Convert to DESeq2 object
diagadds <- phyloseq_to_deseq2(ps_mock_analyze, ~TRT)

# Estimate size factors using geometric means
gm_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans <- apply(counts(diagadds), 1, gm_mean)
diagadds <- estimateSizeFactors(diagadds, geoMeans = geoMeans)

# Run DESeq2
diagadds <- DESeq(diagadds, fitType = "local")

# Extract results
res_DESeq2 <- results(diagadds, cooksCutoff = FALSE)
alpha <- 0.05
sigtab_0.05 <- res_DESeq2[which(res_DESeq2$padj < alpha), ]

# Add taxonomy and ASV ID
sigtab_0.05_COL_Control <- cbind(
  as(sigtab_0.05, "data.frame"),
  as(tax_table(ps_mock_analyze)[rownames(sigtab_0.05), ], "matrix")
)
sigtab_0.05_COL_Control$ASV_ID <- rownames(sigtab_0.05_COL_Control)

# Select top 20 by absolute log2FoldChange
top_sig <- sigtab_0.05_COL_Control[order(abs(sigtab_0.05_COL_Control$log2FoldChange), decreasing = TRUE), ][1:20, ]

# Clean taxonomic labels
top_sig$Genus[is.na(top_sig$Genus)] <- top_sig$Family[is.na(top_sig$Genus)]
top_sig$Species[is.na(top_sig$Species)] <- ""

# Create display name: Genus species [ASV_ID]
top_sig$Taxon <- paste0(top_sig$Genus, " ", top_sig$Species, " [", top_sig$ASV_ID, "]")

# Order by fold change for plot axis
top_sig$Taxon <- factor(top_sig$Taxon, levels = top_sig$Taxon[order(top_sig$log2FoldChange)])

# Add regulation direction
top_sig$Regulation <- ifelse(top_sig$log2FoldChange > 0, "Up in TL", "Down in TL")

# Plot the barplot
ggplot(top_sig, aes(x = Taxon, y = log2FoldChange, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(log2FoldChange, 1)), 
            hjust = ifelse(top_sig$log2FoldChange > 0, -0.1, 1.1), 
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Up in TL" = "#00BFC4", "Down in TL" = "#F8766D")) +
  labs(
    title = "Top 20 Differentially Abundant ASVs (TL vs Control)",
    y = "Log2 Fold Change (TL vs Control)",
    x = "Taxa (Genus Species [ASV ID])",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  expand_limits(y = c(min(top_sig$log2FoldChange) - 2, max(top_sig$log2FoldChange) + 2))
