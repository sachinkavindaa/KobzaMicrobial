# Set working directory to your HCC work folder
setwd("/work/samodha/sachin/Kobza")

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
filt_path <- file.path(fastq_files, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))

# Filtering and trimming
out <- filterAndTrim(fnFs, filtFs, truncLen = 240, maxN = 0, maxEE = 2,
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

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
derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# Denoise
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaFs[1]

# Save intermediate objects
saveRDS(derepFs, "/work/samodha/sachin/Kobza/Save_Model/derepFs_after_trimming.rds")
saveRDS(errF, "/work/samodha/sachin/Kobza/Save_Model/errF_after_trimming.rds")
saveRDS(dadaFs, "/work/samodha/sachin/Kobza/Save_Model/dadaFs_after_trimming.rds")

# Sequence table
seqtab <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtab)))
table(nchar(colnames(seqtab)))

# Filter for 240bp
seqtab <- seqtab[, nchar(colnames(seqtab)) == 240]
dim(seqtab)
sum(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
cat("Non-chimeric retention rate:", round(sum(seqtab.nochim)/sum(seqtab) * 100, 2), "%\n")

# Save non-chimeric table
saveRDS(seqtab.nochim, "/work/samodha/sachin/Kobza/Save_Model/seqtab.nochim_after_trimming.rds")
write.table(t(seqtab.nochim),
            "/work/samodha/sachin/Kobza/Txt_Files/seqtab.nochim_table_after_trimming.txt",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Read tracking
getN <- function(x) sum(getUniques(x))
track_seq <- cbind(
  out,
  sapply(dadaFs, getN),
  rowSums(seqtab.nochim)
)
colnames(track_seq) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track_seq) <- sample.names

write.table(track_seq,
            "/work/samodha/sachin/Kobza/Txt_Files/tracking_reads_through_pipeline_after_trimming.txt",
            sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)



seqtab.nochim <- readRDS("/work/samodha/sachin/Kobza/Save_Model/seqtab.nochim_after_trimming.rds")

# Taxonomic assignment
fastaRef <- "/work/samodha/sachin/Kobza/silva_nr99_v138.2_toSpecies_trainset.fa"
taxTab <- assignTaxonomy(seqtab.nochim, refFasta = fastaRef, multithread = TRUE)
saveRDS(taxTab, "/work/samodha/sachin/Kobza/Save_Model/taxTab_after_trimming.rds")





