library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(decontam); packageVersion("decontam")
library(genefilter); packageVersion("genefilter")
library(Biostrings); packageVersion("Biostrings")
library(DECIPHER)
library(phangorn)
library(pairwiseAdonis)
library(vegan)
library(microbiome)
library(compositions)
library(caret)
library(GUniFrac)
library(phytools)
library(ggplot2)
library(ggtree)
library(treeio)
library(ape)     
library(DESeq2)
library(ashr)
library(scales)
library(ggtext)


load("~/Documents/Kobza/250/ps_mock_analyze_new.rds")
ps_mock_analyze
saveRDS(ps_mock_analyze, file = "~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")
ps_mock_analyze <- readRDS("~/Documents/Kobza/250/ps_mock_analyze_fixed.rds")


sample_data(ps_mock_analyze)$TRT <- factor(sample_data(ps_mock_analyze)$TRT, 
                                           levels = c("Control", "Corn Oil"))
View(ps_mock_analyze@sam_data)


#ps_mock_analyze_con_col <- subset_samples(ps_mock_analyze, TRT == "Control" | TRT == "COL")
#ps_mock_analyze_con_col
#View(ps_mock_analyze_con_col@sam_data)

table(is.na(sample_data(ps_mock_analyze)$TRT))
ps_mock_analyze <- subset_samples(ps_mock_analyze, !is.na(TRT))

diagadds = phyloseq_to_deseq2(ps_mock_analyze, ~TRT)


gm_mean = function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans = apply(counts(diagadds), 1, gm_mean)

diagadds = estimateSizeFactors(diagadds, geoMeans = geoMeans)

diagadds = DESeq(diagadds, fitType = "local")

res_DESeq2 <- results(diagadds, cooksCutoff = FALSE)

alpha = 0.05
sigtab_0.05 <- res_DESeq2[which(res_DESeq2$padj < alpha), ]

sigtab_0.05_COL_Control <- cbind(
  as(sigtab_0.05, "data.frame"),
  as(phyloseq::tax_table(ps_mock_analyze)[rownames(sigtab_0.05), ], "matrix")
)

sigtab_0.05_COL_Control$ASV_ID <- rownames(sigtab_0.05_COL_Control)

View(sigtab_0.05_COL_Control)


top_sig <- sigtab_0.05_COL_Control[order(abs(sigtab_0.05_COL_Control$log2FoldChange), decreasing = TRUE), ][1:20, ]
top_sig$Genus[is.na(top_sig$Genus)] <- top_sig$Family[is.na(top_sig$Genus)]
top_sig$Species[is.na(top_sig$Species)] <- ""

# Create a unique and readable label: Genus species [ASV_ID]
top_sig$Taxon <- paste0(top_sig$Genus, " ", top_sig$Species, " [", top_sig$ASV_ID, "]")

# Reorder factor levels based on log2 fold change
top_sig$Taxon <- factor(top_sig$Taxon, levels = top_sig$Taxon[order(top_sig$log2FoldChange)])

# Add direction for coloring (up/down in COL)
top_sig$Regulation <- ifelse(top_sig$log2FoldChange > 0, "Up in COL", "Down in COL")

# Plot
ggplot(top_sig, aes(x = Taxon, y = log2FoldChange, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(log2FoldChange, 1)), 
            hjust = ifelse(top_sig$log2FoldChange > 0, -0.1, 1.1), 
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Up in COL" = "#00BFC4", "Down in COL" = "#F8766D")) +
  labs(
    title = "Top 20 Differentially Abundant ASVs (Corn Oil vs Control)",
    y = "Log2 Fold Change (Corn Oil vs Control)",
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
  expand_limits(y = c(min(top_sig$log2FoldChange) - 5, max(top_sig$log2FoldChange) + 5))
