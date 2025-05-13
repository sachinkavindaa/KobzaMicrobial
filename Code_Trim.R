library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(decontam); packageVersion("decontam")
library(genefilter); packageVersion("genefilter")
library(Biostrings); packageVersion("Biostrings")

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

# 