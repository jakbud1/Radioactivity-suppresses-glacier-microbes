# DADA2, version 1.26 ("benjjneb/dada2", ref = "v1.26")

### Global options ---------------------------------------------------------
set.seed(123)
library(dada2)
library(dplyr)
library(data.table)

getN <- function(x) sum(getUniques(x))

### Run1 (2022-011) --------------------------------------------------------
## Read files
path_V9_r1 <- "Input/Seq/V9/Run1_011"
files_V9_r1 <- list.files("Input/Seq/V9/Run1_011", pattern = ".fastq", full.names = TRUE)

sample.names_V9_r1 <- substr(files_V9_r1, 23, 35)
 
## Trimming  
filts_V9_r1 <- file.path(path_V9_r1, "filtered", paste0(sample.names_V9_r1, "_filt.fastq.gz"))
names(filts_V9_r1) <- sample.names_V9_r1
 
out_V9_r1 <- filterAndTrim(files_V9_r1, filts_V9_r1,  minLen = 90, maxLen = 180,
                            maxN = 0, maxEE = 1.5, truncQ = 2, rm.phix = FALSE,
                            compress = TRUE, trimLeft = 0, trimRight = 0); head(out_V9_r1, 50)
 
## Train parametric error model  
err_V9_r1 <- learnErrors(filts_V9_r1, randomize = TRUE, multithread = TRUE)
 
png("Output/dada2/V9_r1_error_model.png", width = 800, height = 800)
plotErrors(err_V9_r1, nominalQ = TRUE)
dev.off()
 
## Core sample inference algorithm
dadas_V9_r1 <- dada(filts_V9_r1, err = err_V9_r1, pool = FALSE, multithread = TRUE)
dadas_V9_r1[[1]]
 
seqtab_V9_r1 <- makeSequenceTable(dadas_V9_r1)
dim(seqtab_V9_r1)

table(nchar(getSequences(seqtab_V9_r1)))
 
seqtab.nochim_V9_r1 <- removeBimeraDenovo(seqtab_V9_r1, 
                                           method = "consensus", verbose = TRUE, multithread = TRUE)
dim(seqtab.nochim_V9_r1)
sum(seqtab.nochim_V9_r1)/sum(seqtab_V9_r1)
 
## Track the data
track_V9_r1 <- cbind(out_V9_r1, sapply(dadas_V9_r1, getN), rowSums(seqtab.nochim_V9_r1))
 
colnames(track_V9_r1) <- c("input", "filtered", "denoised", "nonchim")
rownames(track_V9_r1) <- sample.names_V9_r1
head(track_V9_r1, 50)

### Run2 (2022-016) ------------------------------------------------------
## Read files
path_V9_r2 <- "Input/Seq/V9/Run2_016"
files_V9_r2 <- list.files("Input/Seq/V9/Run2_016", pattern = ".fastq", full.names = TRUE)
 
sample.names_V9_r2 <- substr(files_V9_r2, 23, 35)  
 
## Trimming  
filts_V9_r2 <- file.path(path_V9_r2, "filtered", paste0(sample.names_V9_r2, "_filt.fastq.gz"))
names(filts_V9_r2) <- sample.names_V9_r2
 
out_V9_r2 <- filterAndTrim(files_V9_r2, filts_V9_r2,  minLen = 90, maxLen = 180,
                            maxN = 0, maxEE = 1.5, truncQ = 2, rm.phix = FALSE,
                           compress = TRUE, trimLeft = 0, trimRight = 0); head(out_V9_r2, 50)
 
## Train parametric error model  
err_V9_r2 <- learnErrors(filts_V9_r2, randomize = TRUE, multithread = TRUE)
 
png("Output/dada2/V9_r2_error_model.png", width = 800, height = 800)
plotErrors(err_V9_r2, nominalQ = TRUE)
dev.off()
 
## Core sample inference algorithm
dadas_V9_r2 <- dada(filts_V9_r2, err = err_V9_r2, pool = FALSE, multithread = TRUE)
dadas_V9_r2[[1]]
 
seqtab_V9_r2 <- makeSequenceTable(dadas_V9_r2)
dim(seqtab_V9_r2)
 
table(nchar(getSequences(seqtab_V9_r2)))
 
seqtab.nochim_V9_r2 <- removeBimeraDenovo(seqtab_V9_r2, 
                                           method = "consensus", verbose = TRUE, multithread = TRUE)
dim(seqtab.nochim_V9_r2)
sum(seqtab.nochim_V9_r2)/sum(seqtab_V9_r2)
 
## Track the data
track_V9_r2 <- cbind(out_V9_r2, sapply(dadas_V9_r2, getN), rowSums(seqtab.nochim_V9_r2))
 
colnames(track_V9_r2) <- c("input", "filtered", "denoised", "nonchim")
rownames(track_V9_r2) <- sample.names_V9_r2
head(track_V9_r2, 50)

### merging data for both runs, assigning taxonomy ------------------------
seqtab.nochim_V9_all <- mergeSequenceTables(seqtab.nochim_V9_r1, seqtab.nochim_V9_r2)

taxa_V9_all <- assignTaxonomy(seqtab.nochim_V9_all, "Input/Seq/V9/silva_128.18s.99_rep_set.dada2.fa.gz",
                              tryRC = TRUE, minBoot = 60, multithread = TRUE)

colnames(taxa_V9_all) <- make.unique(colnames(taxa_V9_all))

taxa.print_V9_all <- taxa_V9_all # Removing sequence rownames for display only
rownames(taxa.print_V9_all) <- NULL
head(taxa.print_V9_all, 100)

### Save data
seq.abundancesV9_all <- setDT(as.data.frame(as.data.frame(t(seqtab.nochim_V9_all))), keep.rownames = TRUE)

colnames(seq.abundancesV9_all)[colnames(seq.abundancesV9_all) == c("rn")] <- c("ASV_seq")

seq.abundancesV9_all <- setDT(as.data.frame(seq.abundancesV9_all), keep.rownames = TRUE)
seq.abundancesV9_all$rn <- paste0("ASV_", seq.abundancesV9_all$rn)
colnames(seq.abundancesV9_all)[colnames(seq.abundancesV9_all) == c("rn")] <- c("ASV")

taxa1.silva_V9 <- setDT(as.data.frame(taxa_V9_all), keep.rownames = TRUE)
colnames(taxa1.silva_V9)[colnames(taxa1.silva_V9) == c("rn")] <- c("ASV_seq")
taxa.abundance.silva_V9 <- merge(seq.abundancesV9_all, taxa1.silva_V9, by = "ASV_seq")
write.csv(data.table(taxa.abundance.silva_V9), "Output/dada2/taxa.silva_V9_all_boot60.csv", row.names = FALSE)

track_V9_all <- data.frame(cbind(rbind(track_V9_r1, track_V9_r2), c(rep("r1",nrow(track_V9_r1)), rep("r2",nrow(track_V9_r2)))))
track_V9_all <- cbind(rownames(track_V9_all), track_V9_all)
colnames(track_V9_all)[colnames(track_V9_all) == c("V5")] <- "run"

track_V9_all <- data.frame(sample = rownames(track_V9_all), track_V9_all)
write.csv(track_V9_all, file = "Output/dada2/track_reads_V9.csv", row.names = FALSE)

sessionInfo()