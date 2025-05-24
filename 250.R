# Set working directory to your HCC work folder

# Load necessary libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(decontam); packageVersion("decontam")
library(genefilter); packageVersion("genefilter")
library(Biostrings); packageVersion("Biostrings")

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

ps_long_seqs <- readRDS("~/Documents/Kobza/250/ps_mock_filtered_ASV_seqs.rds")
head(taxa_names(ps_long_seqs))
View(tax_table(ps_long_seqs))

ps_asv_ids <- readRDS("~/Documents/Kobza/250/ps_mock_filtered_ASV_IDs.rds")
head(taxa_names(ps_asv_ids))
View(tax_table(ps_asv_ids))

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


plot_richness(ps_rarefy, x = "TRT", measures = c("Observed", "Shannon")) +
  geom_boxplot(aes(group = TRT), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2)

# Ordination (e.g., PCoA)
ordu <- ordinate(norm_mock, method = "PCoA", distance = "bray")
plot_ordination(norm_mock, ordu, color = "TRT")


#PERMANOVA

# Load necessary packages
library(vegan)
library(phyloseq)
library(pairwiseAdonis)
library(ggplot2)

# Step 1: Extract and clean metadata
meta <- as(sample_data(norm_mock), "data.frame")

# Remove samples with NA in TRT (e.g., controls)
meta_clean <- meta[!is.na(meta$TRT), ]

# Ensure TRT and Period are factors
meta_clean$TRT <- factor(meta_clean$TRT)
meta_clean$Period <- factor(meta_clean$Period)

# Step 2: Extract and match OTU table
otu_matrix <- as(otu_table(norm_mock), "matrix")
if (taxa_are_rows(norm_mock)) {
  otu_matrix <- t(otu_matrix)
}

# Match OTU table to cleaned metadata
otu_clean <- otu_matrix[rownames(meta_clean), ]

# Step 3: Compute Bray-Curtis distance
bray_dist <- vegdist(otu_clean, method = "bray")

# Step 4: Run PERMANOVA (Treatment effect)
adonis_trt <- adonis2(bray_dist ~ TRT, data = meta_clean, permutations = 999)
print(adonis_trt)

# Step 5: PERMANOVA with interaction (Treatment * Period)
adonis_interaction <- adonis2(bray_dist ~ TRT + Period + TRT:Period, data = meta_clean, permutations = 999)
print(adonis_interaction)

interaction <- adonis2(bray_dist ~ TRT * Period, data = meta_clean, permutations = 999, strata = meta_clean$Animal)
print(interaction)

# Step 6: Check homogeneity of dispersion (adonis assumption)
dispersion_test <- betadisper(bray_dist, meta_clean$TRT)
permutest_disp <- permutest(dispersion_test)
print(permutest_disp)

plot(dispersion_test)

# Add boxplots to visualize group dispersions
boxplot(dispersion_test)

# Step 7: Pairwise PERMANOVA
pairwise_result <- pairwise.adonis(otu_clean, factors = meta_clean$TRT, 
                                   sim.method = "bray", p.adjust.m = "bonferroni")
print(pairwise_result)



# Load necessary libraries
library(phyloseq)
library(igraph)
library(ggplot2)

# Create a sample-wise network using Jaccard distance
net <- make_network(norm_mock, type = "samples", distance = "bray", max.dist = 0.6)

library(igraph)
vcount(net)  # Number of vertices
ecount(net)  # Number of edges

sample_data(norm_mock)$Period <- as.factor(sample_data(norm_mock)$Period)

plot_network(net, norm_mock, type = "samples", color = "TRT", shape = "Period")




