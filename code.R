library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(decontam); packageVersion("decontam")
library(genefilter); packageVersion("genefilter")
library(Biostrings); packageVersion("Biostrings")


fastq_files <- "./" #specify path to where the raw '.fastq' files are
list.files(fastq_files)

fnFs <- sort(list.files(fastq_files, pattern="_R1_001.fastq.gz")) #forward reads
fnRs <- sort(list.files(fastq_files, pattern="_R2_001.fastq.gz")) #reverse reads

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(fastq_files, fnFs)
fnRs <- file.path(fastq_files, fnRs)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


fnFs[1:3]
fnRs[1:3] 

# Look at base quality plots
plotQualityProfile(fnFs[2]) # plot shows quality of forward reads
plotQualityProfile(fnRs[1:2])#

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

write.table(out, "~/Documents/Kobz/out_file.txt", sep = "\t", col.names = NA, row.names = T, quote = F) 

plotQualityProfile(filtFs[1:2])

# Data Statistics after Trimming
sum(out[,1]) #total reads 
sum(out[,2]) #total reads 
sum(out[,1]) - sum(out[,2]) #reads lost
sum(out[,2])/sum(out[,1]) # percentage data retained 

derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

saveRDS(derepFs, file = "~/Documents/Kobza/Save Model/derepFs.rds")
saveRDS(errF, file = "~/Documents/Kobza/Save Model/errF.rds")
saveRDS(dadaFs, file = "~/Documents/Kobza/Save Model/dadaFs.rds")




