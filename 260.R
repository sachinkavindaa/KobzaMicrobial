setwd("~/Documents/Kobza")

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
plotQualityProfile(fnFs[1:3])


# Create filtered folder
filt_path <- file.path(fastq_files, "filtered") # folder for quality-filtered reads
if(!file_test("-d", filt_path))  dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))

out <- read.table("~/Documents/Kobza/260/out_file_after_trimming_260.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

print(paste("Retention rate:", round(sum(out[,2]) / sum(out[,1]) * 100, 2), "%"))

# Plot quality after filtering
plotQualityProfile(filtFs[1:3])

# Statistics
cat("Total input reads: ", sum(out[,1]), "\n")
cat("Total filtered reads: ", sum(out[,2]), "\n")
cat("Reads lost: ", sum(out[,1]) - sum(out[,2]), "\n")
cat("Percentage retained: ", round(sum(out[,2])/sum(out[,1]) * 100, 2), "%\n")

#This step dereplicates your filtered reads
#It speeds up processing and reduces redundancy. Instead of analyzing every single read, it only processes the unique sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)

#derepFs <- readRDS("~/Documents/Kobza/260/derepFs_after_trimming_260.rds")
names(derepFs) <- sample.names

errF <- readRDS("~/Documents/Kobza/260/errF_after_trimming_260.rds")
plotErrors(errF, nominalQ = TRUE)

dadaFs <- readRDS("~/Documents/Kobza/260/dadaFs_after_trimming_260.rds")
dadaFs[1]

# Sequence table
seqtab <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtab)))
table(nchar(colnames(seqtab)))

# Filter for 240bp
seqtab <- seqtab[, nchar(colnames(seqtab)) == 260]
dim(seqtab)
sum(seqtab)

# Remove chimeras
#during PCR amplification when two different parent sequences are mistakenly joined together.
#They are not real biological organisms, just PCR errors.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
cat("Non-chimeric retention rate:", round(sum(seqtab.nochim)/sum(seqtab) * 100, 2), "%\n")

# Save non-chimeric table
saveRDS(seqtab.nochim, "~/Documents/Kobza/260/seqtab.nochim_after_trimming_260.rds")
write.table(t(seqtab.nochim),
            "~/Documents/Kobza/250/seqtab.nochim_table_after_trimming_250.txt",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

seqtab.nochim <- readRDS("~/Documents/Kobza/260/seqtab.nochim_after_trimming_260.rds")

# Read tracking
getN <- function(x) sum(getUniques(x))
track_seq <- cbind(
  out,
  sapply(dadaFs, getN),
  rowSums(seqtab.nochim)
)
colnames(track_seq) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track_seq) <- sample.names

taxTab <- readRDS("~/Documents/Kobza/260/taxTab_260.rds")

mapping_file <- read.csv("~/Documents/Kobza/Mapping.csv", header = TRUE)

row.names(mapping_file) <- as.character(mapping_file[,1])

#Make phyloseq object
#ASV (amplicon sequence variant) table produced by DADA2 after removing chimeras.
#metadata table (e.g., sample IDs, groupings, time points).
#A matrix of taxonomic assignments
ps_mock <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                    sample_data(mapping_file),
                    tax_table(taxTab))

saveRDS(ps_mock, "~/Documents/Kobza/260/ps_mock_260.rds")

# Fitering out (removing) potential contaminants using 'decontam' package 
#(The more negative controls you have the better this works!)
#Identifying and optionally removing them, you get cleaner, more biologically relevant data.
sample_data(ps_mock)$is.neg <- sample_data(ps_mock)$TYPE == "NEG_CON"
contamdf.prev <- isContaminant(ps_mock, method = "prevalence", neg = "is.neg", threshold = 0.1) # default threshold of 0.1 used
View(contamdf.prev)

hist(contamdf.prev$p, 100, ylim = c(0,400), xlim = c(0,1))
table(contamdf.prev$contaminant) 
#1024 ASVs were not considered contaminants
#2 ASV was flagged as a contaminant based on the threshold

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
#Your primers may also amplify non-bacterial sequences
#They confound bacterial community analysis
remove_kingdom <- c("Eukaryota", "Archaea") 
ps_mock_final_ps <- subset_taxa(ps_mock_filtered, !Kingdom %in% remove_kingdom & (Family != "Mitochondria" | is.na(Family)))
ps_mock_final_ps 

# remove chloroplasts as well
#Chloroplast and mitochondrial 16S rRNA genes are eukaryotic (plant/host) 
#organelles that resemble bacterial genes due to their evolutionary origin.
#They're not microbial
ps_mock_final_ps <- subset_taxa(ps_mock_final_ps, Order != "Chloroplast" | is.na(Order))

#Create DNA sequence set
dna <- Biostrings::DNAStringSet(taxa_names(ps_mock_final_ps))

#Assign names to the sequences
names(dna) <- taxa_names(ps_mock_final_ps)

#Add the DNA sequences to the phyloseq object
ps_mock_final_ps <- merge_phyloseq(ps_mock_final_ps, dna)

#Rename taxa as ASV1, ASV2
taxa_names(ps_mock_final_ps) <- paste0("ASV", seq(ntaxa(ps_mock_final_ps)))


View(t(ps_mock_final_ps@otu_table))

taxa_names(ps_mock_final_ps)
View(tax_table(ps_mock_final_ps))

#filtering out neg & Pos controls
ps_mock_neg <- subset_samples(ps_mock_final_ps, TYPE != "NEG_CON" & Animal_ID != "PCR positive")
ps_mock_neg

#View(t(ps_mock_neg@otu_table))

save(ps_mock_neg, file = "~/Documents/Kobza/260/ps_mock_neg.rds")

#filtering on prevelance and total abundance to remove singletons and spurious ASVs

#prevalence

#It calculates how many samples each taxon appears in
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

#This keeps only ASVs that appear in more than 1 sample
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
save(ps_mock_analyze, file = "~/Documents/Kobza/260/ps_mock_analyze.rds")
load("~/Documents/Kobza/260/ps_mock_analyze.rds")
sum(otu_table(ps_mock_analyze))

library("vegan")

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
save(norm_mock, file= "~/Documents/Kobza/260/norm_mock.rds")

#ps_mock_neg <- subset_samples(ps_mock_neg, Animal != 'PCR positive')
#sample_sums(ps_mock_neg)

#alpha diversity

#simple example of alpha diversity. be sure to check out the above link for other ways to analyze
set.seed(1234)
ps_rarefy <- rarefy_even_depth(ps_mock_neg, 
                               sample.size = min(sample_sums(ps_mock_neg)),
                               rngseed = T, 
                               replace = TRUE, 
                               trimOTUs = TRUE, 
                               verbose = TRUE)
#ps_mock_neg

#View(ps_mock_neg@sam_data)


shannon_rarefy_mock <- estimate_richness(ps_rarefy, split = TRUE, measures = c("Shannon"))
head(shannon_rarefy_mock)
sample_sums(ps_rarefy) 

obser_rarefy_mock <- estimate_richness(ps_rarefy, split = TRUE, measures = c("Observed"))
head(obser_rarefy_mock)

sample_sums(ps_rarefy)


  # Plot alpha diversity
  plot_richness(ps_rarefy, x = "TRT", measures = c("Shannon", "Observed"), color = "TRT") + theme_bw()
  
  plot_richness(ps_rarefy, x = "TRT", measures = c("Observed", "Shannon")) +
    geom_boxplot(aes(group = TRT), alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2)
  
  
  # Ordination (e.g., PCoA)
  ordu <- ordinate(ps_mock_analyze, method = "PCoA", distance = "bray")
  plot_ordination(ps_mock_analyze, ordu, color = "TRT")

#PERMANOVA

# Load necessary packages
library(vegan)
library(phyloseq)
library(pairwiseAdonis)
library(ggplot2)

# Step 1: Extract and clean metadata
meta <- as(sample_data(ps_mock_analyze), "data.frame")

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
