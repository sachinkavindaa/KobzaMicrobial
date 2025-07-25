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



# Input FASTQ files
fastq_files <- "./"
fnFs <- sort(list.files(fastq_files, pattern = "_R1_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])

out <- read.table("~/Documents/Kobza/250/out_file_after_trimming_250.txt", 
                  header = TRUE, sep = "", stringsAsFactors = FALSE)

# Data Statistics after Trimming
sum(out[,1]) #total reads 
sum(out[,2]) #total reads 
sum(out[,1]) - sum(out[,2]) #reads lost
sum(out[,2])/sum(out[,1]) # percentage data retained 

derepFs <- readRDS("~/Documents/Kobza/250/derepFs_after_trimming_250.rds")
names(derepFs) <- sample.names

errF <- readRDS("~/Documents/Kobza/250/errF_after_trimming_250.rds")
plotErrors(errF, nominalQ = TRUE)


dadaFs <- readRDS("~/Documents/Kobza/250/dadaFs_after_trimming_250.rds")
dadaFs[1]


dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

seqtab <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtab)))
table(nchar(colnames(seqtab)))

seqtab <- seqtab[, nchar(colnames(seqtab)) == 250]
dim(seqtab)
sum(seqtab)



seqtab.nochim <- readRDS("~/Documents/Kobza/250/seqtab.nochim_after_trimming_250.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
cat("Non-chimeric retention rate:", round(sum(seqtab.nochim)/sum(seqtab) * 100, 2), "%\n")

View(t(seqtab.nochim))
dim(seqtab.nochim)


# Read tracking
getN <- function(x) sum(getUniques(x))
track_seq <- cbind(
  out,
  sapply(dadaFs, getN),
  rowSums(seqtab.nochim)
)

colnames(track_seq) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track_seq) <- sample.names

taxTab <- readRDS("~/Documents/Kobza/250/taxTab_250.rds")
mapping_file <- read.csv("~/Documents/Kobza/Mapping.csv", header = TRUE)
row.names(mapping_file) <- as.character(mapping_file[,1])
dim(taxTab)

ps_mock <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                    sample_data(mapping_file),
                    tax_table(taxTab))
ps_mock

saveRDS(ps_mock, "~/Documents/Kobza/250/ps_mock_new.rds")

sample_data(ps_mock)$is.neg <- sample_data(ps_mock)$TYPE == "NEG_CON"
contamdf.prev <- isContaminant(ps_mock, method = "prevalence", neg = "is.neg", threshold = 0.1) # default threshold of 0.1 used
View(contamdf.prev)

hist(contamdf.prev$p, 100, ylim = c(0,400), xlim = c(0,1))
table(contamdf.prev$contaminant) 

ps_mock_no_contam <- prune_taxa(!contamdf.prev$contaminant, ps_mock)
ps_mock_no_contam

taxa_names(ps_mock_no_contam)
sum(otu_table(ps_mock_no_contam))/sum(otu_table(ps_mock)) 

#Each ASV count becomes a fraction of the total — i.e., proportion (between 0 and 1)
ps_mock_no_contam.prop  <- transform_sample_counts(ps_mock_no_contam, function (otu) otu/sum(otu)) #Lambda function
ps_mock_no_contam.prop

# set criteria for filter function
flist <- filterfun(kOverA(1,0.0015)) # criteria are at least 0.15% abundance in at least 2 samples
#k = 2: The ASV must appear in at least 2 samples
#A = 0.0015: And in those samples, it must have relative abundance > 0.0015
#Improves statistical power by removing noise

taxa.to.filter <- filter_taxa(ps_mock_no_contam.prop, flist) #create a list of ASVs that meet flist criteria
str(taxa.to.filter)

ps_mock_filtered.prop <- prune_taxa(taxa.to.filter,ps_mock_no_contam.prop)
ps_mock_filtered.prop

taxa_names(ps_mock_filtered.prop)

ps_mock_filtered <- prune_taxa(taxa_names(ps_mock_filtered.prop),ps_mock_no_contam)
ps_mock_filtered

sum(otu_table(ps_mock_filtered))/sum(otu_table(ps_mock_no_contam)) 

remove_kingdom <- c("Eukaryota", "Archaea") 
ps_mock_final_ps <- subset_taxa(ps_mock_filtered, !Kingdom %in% remove_kingdom & (Family != "Mitochondria" | is.na(Family)))
ps_mock_final_ps 

View((ps_mock_final_ps@otu_table))

ps_mock_final_ps <- subset_taxa(ps_mock_final_ps, Order != "Chloroplast" | is.na(Order))

ps_mock_final_ps <- subset_samples(ps_mock_final_ps, TYPE != "NEG_CON" & Animal_ID != "PCR positive")

tree <- read.tree("~/Downloads/mothur/HMA_INF_ASVs.phylip.tre")
ps_tree <- merge_phyloseq(ps_mock_final_ps, tree)


dna <- Biostrings::DNAStringSet(taxa_names(ps_mock_final_ps))
names(dna) <- taxa_names(ps_mock_final_ps)
ps_mock_asv <- merge_phyloseq(ps_mock, dna)
taxa_names(ps_mock_asv) <- paste0("ASV_", seq(ntaxa(ps_mock_asv)))
ps_mock_asv

taxa_names(ps_mock_final_ps)

taxa_names(ps_mock_analyze)

taxa_names(tree)


# Save as FASTA
writeXStringSet(dna, filepath = "~/Downloads/mothur/HMA_INF_ASVs.fa")

ps_mock_final_ps <- merge_phyloseq(ps_mock_final_ps, dna)
taxa_names(ps_mock_final_ps) <- paste0("ASV", seq(ntaxa(ps_mock_final_ps)))

taxa_names(ps_mock_filtered) <- paste0("ASV_", seq(ntaxa(ps_mock_filtered)))
saveRDS(ps_mock_filtered, "~/Documents/Kobza/250/ps_mock_filtered_ASV_IDs_new.rds")



ps_mock_neg <- subset_samples(ps_mock_final_ps, TYPE != "NEG_CON" & Animal_ID != "PCR positive")
ps_mock_neg



save(ps_mock_neg, file = "~/Documents/Kobza/250/ps_mock_neg_new.rds")

prevdf_ps= apply(X = otu_table(ps_mock_neg), 
                 MARGIN = ifelse(taxa_are_rows(ps_mock_neg), yes = 1, no = 2), 
                 FUN = function(x){sum(x > 0)})

prevdf_ps <- data.frame(Prevalence= prevdf_ps, TotalAbundance=taxa_sums(ps_mock_neg))
View(prevdf_ps)

ps_mock_prev <- rownames(prevdf_ps)[prevdf_ps$Prevalence >= 0] #based off of known positive control that was sequenced; for ease of today and a small data set, only set to 1
ps_mock_prev

ps_mock_prev <- prune_taxa(ps_mock_prev, ps_mock_neg)
ps_mock_prev

sum(otu_table(ps_mock_prev))/sum(otu_table(ps_mock_neg))
sum(otu_table(ps_mock_prev))
save(ps_mock_prev, file = "~/Documents/Kobza/250/ps_mock_analyze_250_after_removed_single.rds")

abund_ps= apply(X = otu_table(ps_mock_prev), 
                MARGIN = ifelse(taxa_are_rows(ps_mock_prev), yes = 1, no = 2), 
                FUN = function(x){sum(x > 0)}) 

abund_ps <- data.frame(Prevalence= abund_ps, TotalAbundance=taxa_sums(ps_mock_prev))
View(abund_ps)

abund_taxa <- rownames(abund_ps)

ps_mock_analyze <- prune_taxa(abund_taxa, ps_mock_prev)
ps_mock_analyze

save(ps_mock_analyze, file = "~/Documents/Kobza/250/ps_mock_analyze_new.rds")

otu_mat <- as(otu_table(ps_mock_analyze), "matrix")
if (taxa_are_rows(ps_mock_analyze)) {
  otu_mat <- t(otu_mat)
}

sample_sums <- rowSums(otu_mat)

rarecurve(otu_mat, step = 50, cex = 0.5)
abline(v = min(sample_sums), col = "red", lty = 2)

View((ps_mock_analyze@otu_table))


norm_mock <-  transform_sample_counts(ps_mock_analyze, function(x) x / sum(x) )
save(norm_mock, file= "~/Documents/Kobza/250/norm_mock_new.rds")

set.seed(1234)
ps_rarefy <- rarefy_even_depth(ps_mock_analyze, 
                               sample.size = min(sample_sums(ps_mock_analyze)),
                               rngseed = T, 
                               replace = TRUE, 
                               trimOTUs = TRUE, 
                               verbose = TRUE)

shannon_rarefy_mock <- estimate_richness(ps_rarefy, split = TRUE, measures = c("Shannon"))
head(shannon_rarefy_mock)
sample_sums(ps_rarefy) 

obser_rarefy_mock <- estimate_richness(ps_rarefy, split = TRUE, measures = c("Observed"))
head(obser_rarefy_mock)

plot_richness(ps_rarefy, x = "TRT", measures = c("Observed", "Shannon"), color = "TRT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############ Beta Diversity #############################



# Perform PCoA ordination
ordu <- ordinate(ps_mock_analyze, method = "PCoA", distance = "bray")

# Extract ordination coordinates
ordu_df <- as.data.frame(ordu$vectors[, 1:2])
ordu_df$SampleID <- rownames(ordu_df)

# Extract and prepare metadata
metadata_df <- as(sample_data(ps_mock_analyze), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# Merge ordination coordinates with metadata
plot_df <- merge(ordu_df, metadata_df, by = "SampleID")

# Calculate percentage variance for axes
percent_var <- round(100 * ordu$values$Relative_eig[1:2], 1)

# Plot using ggplot2 with animal labels
ggplot(plot_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = TRT, shape = Period), size = 3) +
  #geom_path(aes(group = Animal), color = "gray40", linetype = "dashed", alpha = 0.6) +
  geom_text(aes(label = Animal), vjust = -1, size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA of Rumen Microbial Communities",
    x = paste0("Axis 1 (", percent_var[1], "%)"),
    y = paste0("Axis 2 (", percent_var[2], "%)"),
    color = "Treatment",
    shape = "Sampling Period"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

dist_mat_clean <- phyloseq::distance(ps_mock_analyze, method = "bray")
permanova <- adonis2(dist_mat_clean ~ Animal + TRT + Period, data = metadata_df, by = "terms")
print(permanova)



########################## PCA (Atchison) #########################

# CLR transformation using the microbiome package
ps_clr <- microbiome::transform(ps_mock_analyze, "clr")

# Extract CLR-transformed OTU matrix
otu_mat <- as(otu_table(ps_clr), "matrix")
if (taxa_are_rows(ps_clr)) {
  otu_mat <- t(otu_mat)
}

# PCA on CLR-transformed data
pca_result <- prcomp(otu_mat, scale. = TRUE)

# Get first 2 PCs and variance explained
pca_df <- as.data.frame(pca_result$x[, 1:2])
percent_var <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Add sample IDs
pca_df$SampleID <- rownames(pca_df)

# Metadata for plotting and PERMANOVA
metadata_df <- as(sample_data(ps_mock_analyze), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# Merge PCA and metadata
plot_df <- merge(pca_df, metadata_df, by = "SampleID")

# Plot PCA with ggplot2
ggplot(plot_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = TRT, shape = Period), size = 3) +
  #geom_path(aes(group = Animal), color = "gray40", linetype = "dashed", alpha = 0.6) +
  geom_text(aes(label = Animal), vjust = -1, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Rumen Microbial Communities (Aitchison's)",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Treatment",
    shape = "Sampling Period"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# PERMANOVA based on Euclidean distance of CLR-PCA input
dist_mat_pca <- dist(otu_mat, method = "euclidean")
permanova_pca <- adonis2(dist_mat_pca ~ Animal + TRT + Period, data = metadata_df, by = "terms")
print(permanova_pca)


############################## Unifrac ########################################


run_unifrac_pcoa <- function(ps_obj, method, title_label) {
  # 1. PCoA ordination
  ordu <- ordinate(ps_obj, method = "PCoA", distance = method)
  
  # 2. Ordination coordinates
  ordu_df <- as.data.frame(ordu$vectors[, 1:2])
  ordu_df$SampleID <- rownames(ordu_df)
  
  # 3. Metadata
  metadata_df <- as(sample_data(ps_obj), "data.frame")
  metadata_df$SampleID <- rownames(metadata_df)
  
  # 4. Merge
  plot_df <- merge(ordu_df, metadata_df, by = "SampleID")
  
  # 5. Percent variance
  percent_var <- round(100 * ordu$values$Relative_eig[1:2], 1)
  
  # 6. Plot
  p <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = TRT, shape = Period), size = 3) +
    #geom_path(aes(group = Animal), color = "gray40", linetype = "dashed", alpha = 0.6) +
    geom_text(aes(label = Animal), vjust = -1, size = 3) +
    theme_minimal() +
    labs(
      title = title_label,
      x = paste0("Axis 1 (", percent_var[1], "%)"),
      y = paste0("Axis 2 (", percent_var[2], "%)"),
      color = "Treatment",
      shape = "Sampling Period"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  
  # 7. PERMANOVA
  dist_mat <- phyloseq::distance(ps_obj, method = method)
  permanova <- adonis2(dist_mat ~ Animal + TRT + Period, data = metadata_df, by = "terms")
  print(permanova)
}



# Unweighted UniFrac
run_unifrac_pcoa(ps_tree, method = "unifrac", title_label = "Unweighted UniFrac PCoA")

# Weighted UniFrac
run_unifrac_pcoa(ps_tree, method = "wunifrac", title_label = "Weighted UniFrac PCoA")




sample_data(ps_mock_analyze)$TRT <- factor(sample_data(ps_mock_analyze)$TRT, 
                                           levels = c("Control", "TL"))
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
    title = "Top 20 Differentially Abundant ASVs (COL vs Control)",
    y = "Log2 Fold Change (COL vs Control)",
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





