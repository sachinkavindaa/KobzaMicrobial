setwd("~/Documents/Kobza")

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(decontam); packageVersion("decontam")
library(genefilter); packageVersion("genefilter")
library(Biostrings); packageVersion("Biostrings")
library("vegan")

fastq_files <- "./" #specify path to where the raw '.fastq' files are
list.files(fastq_files)

fnFs <- sort(list.files(fastq_files, pattern="_R1_001.fastq.gz")) #forward reads

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(fastq_files, fnFs)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

fnFs[1:3]

# Look at base quality plots
plotQualityProfile(fnFs[1:3]) # plot shows quality of forward reads


filt_path <- file.path(fastq_files, "filtered") # folder for quality-filtered reads
if(!file_test("-d", filt_path))  dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))

out <- filterAndTrim(
  fnFs, filtFs, 
  truncLen = c(240), # R1=240, R2=230 (after trimming 10 cycles)
  maxN = 0, 
  maxEE = 2,        # Moderate error thresholds
  truncQ = 2,             # Truncate at Q≤2 (less aggressive than Q=11)
  rm.phix = TRUE,
  compress = TRUE, 
  multithread = TRUE
)

print(paste("Retention rate:", round(sum(out[,2])/sum(out[,1])*100, 2), "%"))

write.table(out, "~/Documents/Kobza/Txt Files/out_file_after_trimming.txt", sep = "\t", col.names = NA, row.names = T, quote = F) 

plotQualityProfile(filtFs[1:3])

# Data Statistics after Trimming
sum(out[,1]) #total reads 
sum(out[,2]) #total reads 
sum(out[,1]) - sum(out[,2]) #reads lost
sum(out[,2])/sum(out[,1]) # percentage data retained 

derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names

errF <- learnErrors(filtFs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs[1] #checking how many sequence variants were inferred from first sample


saveRDS(derepFs, file = "~/Documents/Kobza/Save Model/derepFs_after_trimming.rds")
saveRDS(errF, file = "~/Documents/Kobza/Save Model/errF_after_trimming.rds")
saveRDS(dadaFs, file = "~/Documents/Kobza/Save Model/dadaFs_after_trimming.rds")

# Construct sample-by-sequence table
seqtab <- makeSequenceTable(dadaFs)

# Check sequence length distribution
table(nchar(getSequences(seqtab)))
sum(seqtab)

# Check the sequence lengths
table(nchar(colnames(seqtab)))

# Filter to 240 bp sequences only
seqtab <- seqtab[, nchar(colnames(seqtab)) == 240]

# Check results
dim(seqtab)
sum(seqtab)

#Identification and removal of chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, "~/Documents/Kobza/Save Model/seqtab.nochim_after_trimming.rds")
write.table(t(seqtab.nochim), "~/Documents/Kobza/Txt Files/seqtab.nochim_table_aftter_trimming.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = F)
seqtab.nochim <- readRDS("~/Documents/Kobza/Save Model/seqtab.nochim_after_trimming.rds")
# track reads through the process 
getN <- function(x) sum(getUniques(x))

track_seq <- cbind(
  out,                             # filtered reads
  sapply(dadaFs, getN),            # denoised forward reads
  rowSums(seqtab.nochim)           # non-chimeric reads
)

colnames(track_seq) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track_seq) <- sample.names
head(track_seq)
write.table(track_seq, '~/Documents/Kobza/Txt Files/tracking_reads_through_pipeline_both_sets_after_trimming.txt', sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)

fastaRef <- "silva_nr99_v138.2_toSpecies_trainset.fa" 
taxTab <- assignTaxonomy(seqtab.nochim, refFasta = fastaRef, multithread = TRUE)

taxTab <- readRDS("~/Documents/Kobza/taxTab.rds")

mapping_file <- read.csv("~/Documents/Kobza/Mapping_File.csv", header = TRUE)

row.names(mapping_file) <- as.character(mapping_file[,1])

#Make phyloseq object

ps_mock <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                    sample_data(mapping_file),
                    tax_table(taxTab))

saveRDS(ps_mock, "~/Documents/Kobza/Save Model/ps_mock.rds")

## Fitering out (removing) potential contaminants using 'decontam' package 
#(The more negative controls you have the better this works!)
sample_data(ps_mock)$is.neg <- sample_data(ps_mock)$Type == "NEG_CON"
contamdf.prev <- isContaminant(ps_mock, method = "prevalence", neg = "is.neg", threshold = 0.1) # default threshold of 0.1 used
View(contamdf.prev)

hist(contamdf.prev$p, 100, ylim = c(0,400), xlim = c(0,1))
table(contamdf.prev$contaminant) 

ps_mock_no_contam <- prune_taxa(!contamdf.prev$contaminant, ps_mock)
ps_mock_no_contam 

sum(otu_table(ps_mock_no_contam))/sum(otu_table(ps_mock)) 

## Filtering out low abundant/prevalent ASVs using 
#criteria proposed by Tom et al.(https://www.journalofdairyscience.org/article/S0022-0302(24)00961-5/fulltext) 

# Need to convert to a relative abundance table
ps_mock_no_contam.prop  <- transform_sample_counts(ps_mock_no_contam, function (otu) otu/sum(otu))
ps_mock_no_contam.prop

# set criteria for filter function
flist <- filterfun(kOverA(2,0.0015)) # criteria are at least 0.15% abundance in at least 2 samples

# find out which taxa need to be filtered based on these criteria
taxa.to.filter <- filter_taxa(ps_mock_no_contam.prop, flist) #create a list of ASVs that meet flist criteria
str(taxa.to.filter)

ps_mock_filtered.prop <- prune_taxa(taxa.to.filter,ps_mock_no_contam.prop)
ps_mock_filtered.prop

taxa_names(ps_mock_filtered.prop)

ps_mock_filtered <- prune_taxa(taxa_names(ps_mock_filtered.prop),ps_mock_no_contam)

sum(otu_table(ps_mock_filtered))/sum(otu_table(ps_mock_no_contam)) 

## Removing non-target taxa (Archaea and Eukaryota) 
# filtering out Eukaryota, Archaea, and Mitochondria sequences
remove_kingdom <- c("Eukaryota", "Archaea") 
ps_mock_final_ps <- subset_taxa(ps_mock_filtered, !Kingdom %in% remove_kingdom & (Family != "Mitochondria" | is.na(Family)))
ps_mock_final_ps 

# replace the long ASV sequences with shorter IDs (ASV_1, ASV_2, ASV_3)
saveRDS(ps_mock_filtered, "~/Documents/Kobza/Save Model/ps_mock_filtered_ASV_seqs.rds")

taxa_names(ps_mock_filtered) <- paste0("ASV_", seq(ntaxa(ps_mock_filtered)))

saveRDS(ps_mock_filtered, "~/Documents/Kobza/Save Model/ps_mock_filtered_ASV_IDs.rds")

ps_long_seqs <- readRDS("~/Documents/Kobza/Save Model/ps_mock_filtered_ASV_seqs.rds")
head(taxa_names(ps_long_seqs))
View(tax_table(ps_long_seqs))

ps_asv_ids <- readRDS("~/Documents/Kobza/Save Model/ps_mock_filtered_ASV_IDs.rds")
head(taxa_names(ps_asv_ids))
View(tax_table(ps_asv_ids))

#filtering out neg controls
ps_mock_neg <- subset_samples(ps_mock_final_ps, Type != "NEG_CON")
ps_mock_neg #26 samples with 823 taxa

save(ps_mock_neg, file = "~/Documents/Kobza/Save Model/ps_mock_neg.rds")

#filtering on prevelance and total abundance to remove singletons and spurious ASVs

#prevalence
prevdf_ps= apply(X = otu_table(ps_mock_neg), 
                 MARGIN = ifelse(taxa_are_rows(ps_mock_neg), yes = 1, no = 2), 
                 FUN = function(x){sum(x > 0)})

prevdf_ps <- data.frame(Prevalence= prevdf_ps, TotalAbundance=taxa_sums(ps_mock_neg))
View(prevdf_ps)

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

ps_mock_total_abund <- rownames(abund_ps)[abund_ps$TotalAbundance > 100] ##based off of known positive control that was sequenced; for ease of today and a small data set, only set to 100
ps_mock_total_abund

ps_mock_analyze <- prune_taxa(ps_mock_total_abund, ps_mock_prev)
ps_mock_analyze #90 taxa; 19 samples
sum(otu_table(ps_mock_analyze))/sum(otu_table(ps_mock_prev)) #98% 
save(ps_mock_analyze, file = "~/Documents/Kobza/Save Model/ps_mock_analyze.rds")
load("~/Documents/Kobza/Save Model/ps_mock_analyze.rds")
sum(otu_table(ps_mock_analyze)) #total reads 576631 reads

#rarefraction curve
otu_mat <- as(otu_table(ps_mock_analyze), "matrix")
rarecurve(otu_mat, step = 50, cex = 0.5)

#normalize data on a proportional basis for further analysis (minus alpha diversity)
norm_mock <-  transform_sample_counts(ps_mock_analyze, function(x) x / sum(x) )
save(norm_mock, file= "~/Documents/Kobza/Save Model/norm_mock.rds")

#alpha diversity

#simple example of alpha diversity. be sure to check out the above link for other ways to analyze
set.seed(1234)
ps_rarefy <- rarefy_even_depth(ps_mock_neg, 
                               sample.size = min(sample_sums(ps_mock_neg)),
                               rngseed = T, 
                               replace = TRUE, 
                               trimOTUs = TRUE, 
                               verbose = TRUE)

# Check new sample depths
sample_sums(ps_rarefy)

# Plot alpha diversity
plot_richness(ps_rarefy, x = "TRT", measures = c("Shannon", "Observed"))

# Ordination (e.g., PCoA)
ordu <- ordinate(ps_rarefy, method = "PCoA", distance = "bray")
plot_ordination(ps_rarefy, ordu, color = "TRT")










# 