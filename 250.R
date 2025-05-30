# Set working directory to your HCC work folder
# Load necessary libraries
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



# Input FASTQ files
fastq_files <- "./"
fnFs <- sort(list.files(fastq_files, pattern = "_R1_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])

# Create filtered folder
filt_path <- file.path(fastq_files, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))

# Filtering and trimming
out <- filterAndTrim(fnFs, filtFs, truncLen = 250, maxN = 0, maxEE = 2,
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)


out <- read.table("~/Documents/Kobza/250/out_file_after_trimming_250.txt", 
                  header = TRUE, sep = "", stringsAsFactors = FALSE)

# Retention rate
print(paste("Retention rate:", round(sum(out[,2]) / sum(out[,1]) * 100, 2), "%"))

# Save stats
write.table(out, "/work/samodha/sachin/Kobza/Txt_Files/out_file_after_trimming.txt",
            sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# Plot quality after filtering
plotQualityProfile(filtFs[1:3])

# Statistics
cat("Total input reads: ", sum(out[,1]), "\n")
cat("Total filtered reads: ", sum(out[,2]), "\n")
cat("Reads lost: ", sum(out[,1]) - sum(out[,2]), "\n")
cat("Percentage retained: ", round(sum(out[,2])/sum(out[,1]) * 100, 2), "%\n")

# Dereplication

#This step dereplicates your filtered reads
#It speeds up processing and reduces redundancy. Instead of analyzing every single read, it only processes the unique sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)

derepFs <- readRDS("~/Documents/Kobza/250/derepFs_after_trimming_250.rds")
names(derepFs) <- sample.names

# Learn error rates
#This estimates the error rates in your sequencing data 
#by modeling how often each type of sequencing error occurs (e.g., A→G, T→C) at each base position.
errF <- learnErrors(filtFs, multithread = TRUE)
errF <- readRDS("~/Documents/Kobza/250/errF_after_trimming_250.rds")

plotErrors(errF, nominalQ = TRUE)

#Denoise
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaFs <- readRDS("~/Documents/Kobza/250/dadaFs_after_trimming_250.rds")

dadaFs[1]

# Sequence table
seqtab <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtab)))
table(nchar(colnames(seqtab)))

# Filter for 240bp
seqtab <- seqtab[, nchar(colnames(seqtab)) == 250]
dim(seqtab)
sum(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
cat("Non-chimeric retention rate:", round(sum(seqtab.nochim)/sum(seqtab) * 100, 2), "%\n")

# Save non-chimeric table
saveRDS(seqtab.nochim, "~/Documents/Kobza/250/seqtab.nochim_after_trimming_250.rds")
write.table(t(seqtab.nochim),
            "~/Documents/Kobza/250/seqtab.nochim_table_after_trimming_250.txt",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

seqtab.nochim <- readRDS("~/Documents/Kobza/250/seqtab.nochim_after_trimming_250.rds")

# Read tracking
getN <- function(x) sum(getUniques(x))
track_seq <- cbind(
  out,
  sapply(dadaFs, getN),
  rowSums(seqtab.nochim)
)
colnames(track_seq) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track_seq) <- sample.names


# Taxonomic assignment
fastaRef <- "/work/samodha/sachin/Kobza/silva_nr99_v138.2_toSpecies_trainset.fa"
taxTab <- assignTaxonomy(seqtab.nochim, refFasta = fastaRef, multithread = TRUE)
saveRDS(taxTab, "/work/samodha/sachin/Kobza/Save_Model/taxTab_after_trimming.rds")

taxTab <- readRDS("~/Documents/Kobza/250/taxTab_250.rds")

mapping_file <- read.csv("~/Documents/Kobza/Mapping.csv", header = TRUE)

row.names(mapping_file) <- as.character(mapping_file[,1])


#Make phyloseq object
#ASV (amplicon sequence variant) table produced by DADA2 after removing chimeras.
#metadata table (e.g., sample IDs, groupings, time points).
#A matrix of taxonomic assignments
ps_mock <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                    sample_data(mapping_file),
                    tax_table(taxTab))

saveRDS(ps_mock, "~/Documents/Kobza/250/ps_mock.rds")

# Fitering out (removing) potential contaminants using 'decontam' package 
#(The more negative controls you have the better this works!)
#Identifying and optionally removing them, you get cleaner, more biologically relevant data.
sample_data(ps_mock)$is.neg <- sample_data(ps_mock)$TYPE == "NEG_CON"
contamdf.prev <- isContaminant(ps_mock, method = "prevalence", neg = "is.neg", threshold = 0.1) # default threshold of 0.1 used
View(contamdf.prev)

hist(contamdf.prev$p, 100, ylim = c(0,400), xlim = c(0,1))
table(contamdf.prev$contaminant) 
#817 ASVs were not considered contaminants
#1 ASV was flagged as a contaminant based on the threshold

#Remove the contaminanted 
ps_mock_no_contam <- prune_taxa(!contamdf.prev$contaminant, ps_mock)
ps_mock_no_contam 

sum(otu_table(ps_mock_no_contam))/sum(otu_table(ps_mock)) 

## Filtering out low abundant/prevalent ASVs using 
#criteria proposed by Tom et al.(https://www.journalofdairyscience.org/article/S0022-0302(24)00961-5/fulltext) 

#Need to convert to a relative abundance table
#Normalizing data to relative abundance so that each sample is scaled equally
#Each ASV count becomes a fraction of the total — i.e., proportion (between 0 and 1)
ps_mock_no_contam.prop  <- transform_sample_counts(ps_mock_no_contam, function (otu) otu/sum(otu)) #Lambda function
ps_mock_no_contam.prop

# set criteria for filter function
flist <- filterfun(kOverA(2,0.0015)) # criteria are at least 0.15% abundance in at least 2 samples
#k = 2: The ASV must appear in at least 2 samples
#A = 0.0015: And in those samples, it must have relative abundance > 0.0015
#Improves statistical power by removing noise

# find out which taxa need to be filtered based on these criteria
taxa.to.filter <- filter_taxa(ps_mock_no_contam.prop, flist) #create a list of ASVs that meet flist criteria
str(taxa.to.filter)

#Contaminants removed
#Rare, low-abundance ASVs filtered out
ps_mock_filtered.prop <- prune_taxa(taxa.to.filter,ps_mock_no_contam.prop)
ps_mock_filtered.prop

taxa_names(ps_mock_filtered.prop)

#Synchronize your filtered ASVs back to your original count-based dataset.
ps_mock_filtered <- prune_taxa(taxa_names(ps_mock_filtered.prop),ps_mock_no_contam)

sum(otu_table(ps_mock_filtered))/sum(otu_table(ps_mock_no_contam)) 

## Removing non-target taxa (Archaea and Eukaryota) 

# filtering out Eukaryota, Archaea, and Mitochondria sequences
remove_kingdom <- c("Eukaryota", "Archaea") 
ps_mock_final_ps <- subset_taxa(ps_mock_filtered, !Kingdom %in% remove_kingdom & (Family != "Mitochondria" | is.na(Family)))
ps_mock_final_ps 

# remove chloroplasts as well
ps_mock_final_ps <- subset_taxa(ps_mock_final_ps, Order != "Chloroplast" | is.na(Order))

dna <- Biostrings::DNAStringSet(taxa_names(ps_mock_final_ps))
names(dna) <- taxa_names(ps_mock_final_ps)

ps_mock_final_ps <- merge_phyloseq(ps_mock_final_ps, dna)
taxa_names(ps_mock_final_ps) <- paste0("ASV", seq(ntaxa(ps_mock_final_ps)))

saveRDS(ps_mock_filtered, "~/Documents/Kobza/250/ps_mock_filtered_ASV_seqs.rds")

taxa_names(ps_mock_filtered) <- paste0("ASV_", seq(ntaxa(ps_mock_filtered)))

saveRDS(ps_mock_filtered, "~/Documents/Kobza/250/ps_mock_filtered_ASV_IDs.rds")


#filtering out neg controls
ps_mock_neg <- subset_samples(ps_mock_final_ps, TYPE != "NEG_CON" & Animal_ID != "PCR positive")
ps_mock_neg

save(ps_mock_neg, file = "~/Documents/Kobza/250/ps_mock_neg.rds")

#filtering on prevelance and total abundance to remove singletons and spurious ASVs

#prevalence
prevdf_ps= apply(X = otu_table(ps_mock_neg), 
                 MARGIN = ifelse(taxa_are_rows(ps_mock_neg), yes = 1, no = 2), 
                 FUN = function(x){sum(x > 0)})

#Remove ASVs that are only seen once (singletons)
#Remove ASVs that have very low total counts
#Helps you decide which ASVs are rare (low prevalence) or low-abundance, and may be filtered out
prevdf_ps <- data.frame(Prevalence= prevdf_ps, TotalAbundance=taxa_sums(ps_mock_neg))
View(prevdf_ps)


# 1. Get prevalence vector
prev_vals <- apply(
  X = otu_table(ps_mock_neg),
  MARGIN = ifelse(taxa_are_rows(ps_mock_neg), 1, 2),
  FUN = function(x) sum(x > 0)
)

# 2. Combine with total abundance into a proper data frame
prevdf_ps <- data.frame(
  Prevalence = prev_vals,
  TotalAbundance = taxa_sums(ps_mock_neg)
)

# 3. Add ASV names as rownames (in case they are not already)
rownames(prevdf_ps) <- taxa_names(ps_mock_neg)

# 4. Filter ASVs with prevalence > 1
ps_mock_prev <- rownames(prevdf_ps)[prevdf_ps$Prevalence > 1] #based off of known positive control that was sequenced; for ease of today and a small data set, only set to 1
ps_mock_prev

ps_mock_prev <- prune_taxa(ps_mock_prev, ps_mock_neg)
ps_mock_prev

sum(otu_table(ps_mock_prev))/sum(otu_table(ps_mock_neg))

#total abundance  
abund_ps= apply(X = otu_table(ps_mock_prev), 
                MARGIN = ifelse(taxa_are_rows(ps_mock_prev), yes = 1, no = 2), 
                FUN = function(x){sum(x > 0)})

abund_ps <- data.frame(Prevalence= abund_ps, TotalAbundance=taxa_sums(ps_mock_prev))
View(abund_ps)

#Keep only those ASVs where total abundance is > 100
ps_mock_total_abund <- rownames(abund_ps)[abund_ps$TotalAbundance > 100] ##based off of known positive control that was sequenced; for ease of today and a small data set, only set to 100
ps_mock_total_abund

ps_mock_analyze <- prune_taxa(ps_mock_total_abund, ps_mock_prev)
ps_mock_analyze

sum(otu_table(ps_mock_analyze))/sum(otu_table(ps_mock_prev)) #98% 
save(ps_mock_analyze, file = "~/Documents/Kobza/250/ps_mock_analyze.rds")
load("~/Documents/Kobza/250/ps_mock_analyze.rds")
sum(otu_table(ps_mock_analyze))

# Convert to matrix and transpose if taxa are rows
otu_mat <- as(otu_table(ps_mock_analyze), "matrix")
if (taxa_are_rows(ps_mock_analyze)) {
  otu_mat <- t(otu_mat)
}

sample_sums <- rowSums(otu_mat)

# Now plot rarefaction curves
rarecurve(otu_mat, step = 50, cex = 0.5)
abline(v = min(sample_sums), col = "red", lty = 2)

View((ps_mock_analyze@otu_table))

#normalize data on a proportional basis for further analysis (minus alpha diversity)
norm_mock <-  transform_sample_counts(ps_mock_analyze, function(x) x / sum(x) )
save(norm_mock, file= "~/Documents/Kobza/250/norm_mock.rds")

#alpha diversity

#simple example of alpha diversity. be sure to check out the above link for other ways to analyze
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


# Ordination (e.g., PCoA)
ordu <- ordinate(ps_mock_analyze, method = "PCoA", distance = "bray")
plot_ordination(ps_mock_analyze, ordu, color = "TRT", shape = "Period")

percent_var <- round(100 * ordu$values$Relative_eig[1:2], 1)  # PCoA Axis 1 and 2

# Update plot with percentage labels
p <- plot_ordination(ps_mock_analyze, ordu, color = "TRT", shape = "Period") +
  geom_point(size = 3) +
  ggtitle("PCoA of Rumen Microbial Communities") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(
    x = paste0("Axis 1 (", percent_var[1], "%)"),
    y = paste0("Axis 2 (", percent_var[2], "%)"),
    color = "Treatment",
    shape = "Sampling Period"
  )

# Display the plot
print(p)

metadata_df <- as(sample_data(ps_mock_analyze), "data.frame")
dist_mat_clean <- phyloseq::distance(ps_mock_analyze, method = "bray")
permanova <- adonis2(dist_mat_clean ~ TRT + Period, data = metadata_df, by = "terms")
print(permanova)

View(norm_mock@sam_data)








################## Phlogenitic tree #######################

dna_seqs <- refseq(ps_mock_analyze)
alignment <- AlignSeqs(dna_seqs, anchor=NA)

phang_align <- phyDat(as(alignment, "matrix"), type="DNA") # Convert alignment to phangorn format

dm <- dist.ml(phang_align) # Distance matrix
treeNJ <- NJ(dm) # Build initial NJ tree

# Maximum likelihood optimization
fit <- pml(treeNJ, data=phang_align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

tree <- fitGTR$tree

save(tree, file = "~/Documents/Kobza/250/tree.rds")

ps_mock_analyze <- merge_phyloseq(ps_mock_analyze, tree)
norm_tree_mock <- transform_sample_counts(ps_mock_analyze, function(x) x / sum(x))


#################### Prepare and clean metadata #############################

metadata <- as(sample_data(norm_tree_mock), "data.frame")
metadata_clean <- metadata[complete.cases(metadata[, c("TRT", "Period")]), ]
metadata_clean$TRT <- as.factor(metadata_clean$TRT)
metadata_clean$Period <- as.factor(metadata_clean$Period)
samples_to_keep <- rownames(metadata_clean)




run_unifrac_analysis <- function(physeq_obj, method = "unifrac",
                                 metadata_df, sample_ids,
                                 title_label = "PCoA Plot",
                                 save_path = NULL) {
  # 1. Subset and calculate distance matrix
  dist_mat <- phyloseq::distance(physeq_obj, method = method)
  dist_mat_clean <- as.matrix(dist_mat)[sample_ids, sample_ids]
  dist_mat_clean <- as.dist(dist_mat_clean)
  
  # 2. Run PCoA ordination
  ord <- ordinate(physeq_obj, method = "PCoA", distance = method)
  
  # 3. Get % variance explained for axes
  pvar <- round(100 * ord$values$Relative_eig[1:2], 1)
  x_lab <- paste0("PCoA Axis 1 [", pvar[1], "%]")
  y_lab <- paste0("PCoA Axis 2 [", pvar[2], "%]")
  
  # 4. Create plot
  p <- plot_ordination(physeq_obj, ord, color = "TRT", shape = "Period") +
    geom_point(size = 3, alpha = 0.9) +
    ggtitle(title_label) +
    labs(x = x_lab, y = y_lab, color = "Treatment", shape = "Sampling Period") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  
  print(p)

  # 7. Run PERMANOVA
  permanova <- adonis2(dist_mat_clean ~ TRT + Period, data = metadata_df, by = "terms")
  print(permanova)
  
  # 8. Return results
  return(list(plot = p, result = permanova))
}


# Run Unweighted UniFrac PCoA + PERMANOVA
result_unifrac <- run_unifrac_analysis(
  physeq_obj = norm_tree_mock,
  method = "unifrac",
  metadata_df = metadata_clean,
  sample_ids = samples_to_keep,
  title_label = "Unweighted UniFrac",
  save_path = "unweighted_unifrac_pcoa.pdf"  # optional: save as PDF
)

# Run Weighted UniFrac PCoA + PERMANOVA
result_wunifrac <- run_unifrac_analysis(
  physeq_obj = norm_tree_mock,
  method = "wunifrac",
  metadata_df = metadata_clean,
  sample_ids = samples_to_keep,
  title_label = "Weighted UniFrac",
  save_path = "weighted_unifrac_pcoa.pdf"  # optional
)



run_pca_analysis <- function(physeq_obj, metadata_df, sample_ids,
                             title_label = "PCA Ordination",
                             save_path = NULL) {
  # 1. Extract OTU/ASV matrix (samples as rows)
  otu_mat <- as(otu_table(physeq_obj), "matrix")
  if (taxa_are_rows(physeq_obj)) {
    otu_mat <- t(otu_mat)
  }
  
  # 2. Subset to cleaned samples
  otu_mat <- otu_mat[sample_ids, ]
  
  # 3. Perform PCA on scaled data
  pca_result <- prcomp(otu_mat, scale. = TRUE)
  
  # 4. Prepare PCA scores
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$SampleID <- rownames(pca_df)
  
  # 5. Merge with metadata
  meta_df <- metadata_df
  meta_df$SampleID <- rownames(meta_df)
  merged_df <- merge(pca_df, meta_df, by = "SampleID")
  
  # 6. Calculate % variance explained
  pvar <- round(100 * summary(pca_result)$importance[2, 1:2], 1)
  x_lab <- paste0("PC1 (", pvar[1], "%)")
  y_lab <- paste0("PC2 (", pvar[2], "%)")
  
  # 7. Create PCA plot
  p <- ggplot(merged_df, aes(x = PC1, y = PC2, color = TRT, shape = Period)) +
    geom_point(size = 3, alpha = 0.9) +
    ggtitle(title_label) +
    labs(x = x_lab, y = y_lab, color = "Treatment", shape = "Sampling Period") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  print(p)
  
  # 9. Run PERMANOVA on Euclidean distances
  dist_mat <- dist(otu_mat)
  permanova <- adonis2(dist_mat ~ TRT + Period, data = metadata_df, by = "terms")
  print(permanova)
  
  # 10. Return plot and stats
  return(list(plot = p, result = permanova))
}

result_pca <- run_pca_analysis(
  physeq_obj = norm_tree_mock,
  metadata_df = metadata_clean,
  sample_ids = samples_to_keep,
  title_label = "PCA of Rumen Microbial Communities"  # Optional
)



run_aitchison_pca_analysis <- function(physeq_obj,
                                       metadata_df,
                                       sample_ids,
                                       title_label = "Aitchison PCA (CLR-transformed)",
                                       plot_filename = NULL,
                                       verbose = TRUE) {
  # 1. CLR transform
  ps_clr <- microbiome::transform(physeq_obj, "clr")
  
  # 2. Extract OTU/ASV matrix
  otu_mat <- as(otu_table(ps_clr), "matrix")
  if (taxa_are_rows(ps_clr)) {
    otu_mat <- t(otu_mat)
  }
  
  # 3. Subset to valid samples
  otu_mat <- otu_mat[sample_ids, ]
  
  # 4. Run PCA on CLR-transformed data
  pca_result <- prcomp(otu_mat, scale. = TRUE)
  
  # 5. Prepare metadata
  meta_df <- metadata_df
  meta_df$SampleID <- rownames(meta_df)
  
  # 6. Combine PCA scores + metadata
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$SampleID <- rownames(pca_df)
  merged_df <- merge(pca_df, meta_df, by = "SampleID")
  
  # 7. Calculate % variance explained
  pvar <- round(100 * summary(pca_result)$importance[2, 1:2], 1)
  x_lab <- paste0("PC1 (", pvar[1], "%)")
  y_lab <- paste0("PC2 (", pvar[2], "%)")
  
  # 8. Plot
  p <- ggplot(merged_df, aes(x = PC1, y = PC2, color = TRT, shape = Period)) +
    geom_point(size = 3, alpha = 0.9) +
    ggtitle(title_label) +
    labs(x = x_lab, y = y_lab, color = "Treatment", shape = "Sampling Period") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  # 9. Print or save plot
  if (verbose) print(p)
  if (!is.null(plot_filename)) {
    ggsave(plot_filename, p, width = 8, height = 6, dpi = 300)
  }
  
  # 10. Run PERMANOVA (Euclidean on CLR)
  dist_mat <- dist(otu_mat)
  permanova <- adonis2(dist_mat ~ TRT + Period, data = metadata_df, by = "terms")
  if (verbose) print(permanova)
  
  # 11. Return
  return(list(plot = p, result = permanova))
}

result_aitchison <- run_aitchison_pca_analysis(
  physeq_obj = norm_tree_mock,
  metadata_df = metadata_clean,
  sample_ids = samples_to_keep,
  title_label = "Aitchison PCA (CLR-transformed)",
  plot_filename = "figures/aitchison_pca.pdf"  # Save as PDF
)




# Step 1: Extract OTU matrix
otu_mat <- as(otu_table(norm_tree_mock), "matrix")
if (taxa_are_rows(norm_tree_mock)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- otu_mat[samples_to_keep, ]

# Step 2: Convert to compositional object
otu_acomp <- acomp(otu_mat)

# Step 3: ILR transformation
otu_ilr <- ilr(otu_acomp)
otu_ilr_matrix <- as.matrix(otu_ilr)

# Step 4: PCA
pca_result <- prcomp(otu_ilr_matrix, scale. = TRUE)

# Step 4.5: PERMANOVA on ILR Euclidean distance
metadata_df <- metadata_clean
metadata_df$SampleID <- rownames(metadata_df)
ilr_dist <- dist(otu_ilr_matrix)
adonis_ilr <- adonis2(ilr_dist ~ TRT + Period, data = metadata_df, by = "terms")
print(adonis_ilr)

# Step 5: Prepare PCA plot data
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df <- merge(pca_df, metadata_df, by = "SampleID")

# Step 6: Variance explained
pvar <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Step 7: Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = TRT, shape = Period)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA using ILR Transformation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = paste0("PC1 (", pvar[1], "%)"),
    y = paste0("PC2 (", pvar[2], "%)")
  )




sample_ids <- rownames(metadata_clean)
physeq_obj <- norm_tree_mock
physeq_sub <- prune_samples(sample_ids, physeq_obj)

# Root the tree (recommended)
phy_tree(physeq_sub) <- ape::root(phy_tree(physeq_sub), outgroup = taxa_names(physeq_sub)[1], resolve.root = TRUE)

# Compute distance and ordination
dist_mat <- phyloseq::distance(physeq_sub, method = "unifrac")
ord <- ordinate(physeq_sub, method = "PCoA", distance = "unifrac")

group <- metadata_df$TRT[sample_ids]
bd <- betadisper(dist_mat, group)
anova(bd)  # If p < 0.05, dispersion differs across groups



# Ensure consistent metadata and distance matrix
dist_mat <- phyloseq::distance(physeq_sub, method = "unifrac")
meta_sub <- metadata_df[rownames(as.matrix(dist_mat)), ]
group <- meta_sub$TRT

# Run betadisper properly
bd <- betadisper(dist_mat, group)
anova(bd)  # p < 0.05 means dispersion differs significantly across groups


phyloseq::distance(physeq_sub, method = "wunifrac")

# CLR transform
otu_clr <- compositions::clr(t(otu_table(physeq_sub) + 1e-6))  # transpose BEFORE clr

dist_clr <- dist(otu_clr)  # samples in rows
meta_sub <- metadata_clean[rownames(otu_clr), ]  # subset metadata to same samples

# Step 1: CLR on transposed OTU table
otu_clr <- compositions::clr(t(otu_table(physeq_sub) + 1e-6))

# Step 2: Create matching metadata
meta_sub <- metadata_clean[rownames(otu_clr), ]

# Step 3: Remove samples with missing TRT or Period
meta_sub <- meta_sub[complete.cases(meta_sub[, c("TRT", "Period")]), ]

# Step 4: Subset CLR matrix to only valid samples
otu_clr <- otu_clr[rownames(meta_sub), ]

# Step 5: Compute distance and run PERMANOVA
dist_clr <- dist(otu_clr)
adonis2(dist_clr ~ TRT + Period, data = meta_sub)




tree_rooted <- phytools::midpoint.root(phy_tree(norm_tree_mock))

# Replace the tree in your phyloseq object
phy_tree(norm_tree_mock) <- tree_rooted



run_generalized_unifrac <- function(physeq_obj, metadata_df, sample_ids,
                                    alpha_val = 0.5,
                                    title_label = "Generalized UniFrac (d=0.5) PCoA") {
  # Ensure phyloseq is loaded
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Load phyloseq package")
  if (!requireNamespace("GUniFrac", quietly = TRUE)) stop("Load GUniFrac package")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Load ggplot2 package")
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Load vegan package")
  if (!requireNamespace("ape", quietly = TRUE)) stop("Load ape package")
  
  # 1. Extract OTU table
  otu_mat <- as(otu_table(physeq_obj), "matrix")
  if (taxa_are_rows(physeq_obj)) {
    otu_mat <- t(otu_mat)
  }
  otu_mat <- otu_mat[sample_ids, ]
  
  # 2. Get tree
  tree <- phy_tree(physeq_obj)
  
  # Compute Generalized UniFrac distances
  unifracs <- GUniFrac(otu.tab = otu_mat, tree = tree, alpha = c(alpha_val))$unifracs
  dist_mat <- unifracs[, , paste0("d_", alpha_val)]
  dist_mat <- as.dist(dist_mat)
  
  
  # 4. Run PCoA
  pcoa <- ape::pcoa(dist_mat)
  pvar <- round(pcoa$values$Relative_eig[1:2] * 100, 1)
  x_lab <- paste0("PCoA Axis 1 (", pvar[1], "%)")
  y_lab <- paste0("PCoA Axis 2 (", pvar[2], "%)")
  
  # 5. Prepare metadata + ordination data
  metadata_df$SampleID <- rownames(metadata_df)
  ord_df <- as.data.frame(pcoa$vectors[, 1:2])
  ord_df$SampleID <- rownames(ord_df)
  merged_df <- merge(ord_df, metadata_df, by = "SampleID")
  
  # 6. Plot
  p <- ggplot(merged_df, aes(x = Axis.1, y = Axis.2, color = TRT, shape = Period)) +
    geom_point(size = 3, alpha = 0.9) +
    ggtitle(title_label) +
    labs(x = x_lab, y = y_lab, color = "Treatment", shape = "Sampling Period") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  print(p)
  
  # 7. Run PERMANOVA
  permanova <- adonis2(dist_mat ~ TRT + Period, data = metadata_df, by = "terms")
  print(permanova)
  
  # 8. Return
  return(list(plot = p, result = permanova, distance = dist_mat))
}





# 1. Extract OTU table
otu_mat <- as(otu_table(norm_tree_mock), "matrix")
if (taxa_are_rows(norm_tree_mock)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- otu_mat[samples_to_keep, ]

# 2. Get tree
tree <- phy_tree(norm_tree_mock)

# 3. Root the tree if not already
if (!is.rooted(tree)) {
  library(phytools)
  tree <- midpoint.root(tree)
  phy_tree(norm_tree_mock) <- tree
}

# 4. Calculate Generalized UniFrac distance (alpha = 0.5)
unifracs <- GUniFrac(otu.tab = otu_mat, tree = tree, alpha = c(0.5))$unifracs
dist_mat <- unifracs[, , "d_0.5"]
dist_mat <- as.dist(dist_mat)

# 5. Perform PCoA
pcoa <- ape::pcoa(dist_mat)

# 6. Prepare ordination + metadata dataframe
metadata_df <- metadata_clean
metadata_df$SampleID <- rownames(metadata_df)
pcoa_df <- as.data.frame(pcoa$vectors[, 1:2])
pcoa_df$SampleID <- rownames(pcoa_df)
merged_df <- merge(pcoa_df, metadata_df, by = "SampleID")

# 7. Plot with ellipses and centroids
library(ggplot2)
ggplot(merged_df, aes(x = Axis.1, y = Axis.2, color = TRT, shape = Period)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = TRT), linetype = "dashed", level = 0.68, size = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 4, stroke = 2, aes(group = TRT)) +
  labs(
    title = "Generalized UniFrac d(0.5): Rumen Microbiome",
    x = paste0("PCoA Axis 1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA Axis 2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)")
  ) +
  theme_minimal(base_size = 14)


############### Visualise tree ############################
tree <- phy_tree(ps_mock_analyze)

# Get taxonomy table
tax_tab <- as.data.frame(tax_table(ps_mock_analyze))

# Extract ASV tips in the same order as the tree
asv_names <- tree$tip.label

# Choose a taxonomic rank to color (e.g., Phylum)
tax_coloring_rank <- "Phylum"
tax_groups <- tax_tab[asv_names, tax_coloring_rank]

# Replace NAs or empty strings with "Unclassified"
tax_groups[is.na(tax_groups) | tax_groups == ""] <- "Unclassified"

# Assign colors to each group
tax_levels <- unique(tax_groups)
colors <- setNames(brewer.pal(min(length(tax_levels), 8), "Set1"), tax_levels)

# Handle more than 8 groups
if (length(tax_levels) > 8) {
  library(randomcoloR)
  colors <- setNames(distinctColorPalette(length(tax_levels)), tax_levels)
}

tip_colors <- colors[tax_groups]

plot(tree,
     tip.color = tip_colors,
     cex = 0.4,
     label.offset = 0.01,
     main = "Phylogenetic Tree of ASVs by Phylum")

add.scale.bar(length = 0.1, lwd = 1)

# Optional legend
legend("topright", legend = names(colors), col = colors, pch = 19, cex = 0.5, box.lty = 0)
