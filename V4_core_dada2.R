# DADA2, version 1.26 ("benjjneb/dada2", ref = "v1.26")

### Global options ---------------------------------------------------------
set.seed(123)
library(dada2)
library(dplyr)
library(data.table)

getN <- function(x) sum(getUniques(x))

### Run1 (2022-011) --------------------------------------------------------
## Read files
path_V4_r1 <- "Input/Seq/V4/Run1_011"
files_V4_r1 <- list.files("Input/Seq/V4/Run1_011", pattern = ".fastq", full.names = TRUE)

sample.names_V4_r1 <- substr(files_V4_r1, 23, 35)

## Trimming  
filts_V4_r1 <- file.path(path_V4_r1, "filtered", paste0(sample.names_V4_r1, "_filt.fastq.gz"))
names(filts_V4_r1) <- sample.names_V4_r1

out_V4_r1 <- filterAndTrim(files_V4_r1, filts_V4_r1,  minLen = 160, maxLen = 300,
                           maxN = 0, maxEE = 3, truncQ = 2, rm.phix = TRUE,
                           compress = TRUE, trimLeft = 10, trimRight = 25); head(out_V4_r1, 50)

## Train parametric error model  
err_V4_r1 <- learnErrors(filts_V4_r1, randomize = TRUE, multithread = TRUE)

png("Output/dada2/V4_r1_error_model.png", width = 800, height = 800)
plotErrors(err_V4_r1, nominalQ = TRUE)
dev.off()

## Core sample inference algorithm
dadas_V4_r1 <- dada(filts_V4_r1, err = err_V4_r1, pool = FALSE, multithread = TRUE)
dadas_V4_r1[[2]]

seqtab_V4_r1 <- makeSequenceTable(dadas_V4_r1)
dim(seqtab_V4_r1)

table(nchar(getSequences(seqtab_V4_r1)))

seqtab.nochim_V4_r1 <- removeBimeraDenovo(seqtab_V4_r1, 
                                          method = "consensus", verbose = TRUE)
dim(seqtab.nochim_V4_r1)
sum(seqtab.nochim_V4_r1)/sum(seqtab_V4_r1)

## Track the data
track_V4_r1 <- cbind(out_V4_r1, sapply(dadas_V4_r1, getN), rowSums(seqtab.nochim_V4_r1))

colnames(track_V4_r1) <- c("input", "filtered", "denoised", "nonchim")
rownames(track_V4_r1) <- sample.names_V4_r1
head(track_V4_r1, 50)

### Run2 (2022-016) ------------------------------------------------------
## Read files
path_V4_r2 <- "Input/Seq/V4/Run2_016"
files_V4_r2 <- list.files("Input/Seq/V4/Run2_016", pattern = ".fastq", full.names = TRUE)

sample.names_V4_r2 <- substr(files_V4_r2, 23, 35)  

## Trimming  
filts_V4_r2 <- file.path(path_V4_r2, "filtered", paste0(sample.names_V4_r2, "_filt.fastq.gz"))
names(filts_V4_r2) <- sample.names_V4_r2

out_V4_r2 <- filterAndTrim(files_V4_r2, filts_V4_r2,  minLen = 160, maxLen = 300,
                           maxN = 0, maxEE = 3, truncQ = 2, rm.phix = TRUE,
                           compress = TRUE, trimLeft = 10, trimRight = 25); head(out_V4_r2, 50)

## Train parametric error model  
err_V4_r2 <- learnErrors(filts_V4_r2, randomize = TRUE, multithread = TRUE)

png("Output/dada2/V4_r2_error_model.png", width = 800, height = 800)
plotErrors(err_V4_r2, nominalQ = TRUE)
dev.off()

## Core sample inference algorithm
dadas_V4_r2 <- dada(filts_V4_r2, err = err_V4_r2, pool = FALSE, multithread = TRUE)
dadas_V4_r2[[1]]

seqtab_V4_r2 <- makeSequenceTable(dadas_V4_r2)
dim(seqtab_V4_r2)

table(nchar(getSequences(seqtab_V4_r2)))

seqtab.nochim_V4_r2 <- removeBimeraDenovo(seqtab_V4_r2, 
                                          method = "consensus", verbose = TRUE)
dim(seqtab.nochim_V4_r2)
sum(seqtab.nochim_V4_r2)/sum(seqtab_V4_r2)

## Track the data
track_V4_r2 <- cbind(out_V4_r2, sapply(dadas_V4_r2, getN), rowSums(seqtab.nochim_V4_r2))

colnames(track_V4_r2) <- c("input", "filtered", "denoised", "nonchim")
rownames(track_V4_r2) <- sample.names_V4_r2
head(track_V4_r2, 50)

### merging data for both runs, assigning taxonomy -------------------------
seqtab.nochim_V4_all <- mergeSequenceTables(seqtab.nochim_V4_r1, seqtab.nochim_V4_r2)

taxa_V4_all <- assignTaxonomy(seqtab.nochim_V4_all, "Input/Seq/V4/silva_nr99_v138.1_train_set.fa.gz",
                              tryRC = TRUE, minBoot = 95)

taxa.print_V4_all <- taxa_V4_all # Removing sequence rownames for display only
rownames(taxa.print_V4_all) <- NULL
head(taxa.print_V4_all, 100)

### Save data
seq.abundancesV4_all <- setDT(as.data.frame(as.data.frame(t(seqtab.nochim_V4_all))), keep.rownames = TRUE)

colnames(seq.abundancesV4_all)[colnames(seq.abundancesV4_all) == c("rn")] <- c("ASV_seq")

seq.abundancesV4_all <- setDT(as.data.frame(seq.abundancesV4_all), keep.rownames = TRUE)
seq.abundancesV4_all$rn <- paste0("ASV_", seq.abundancesV4_all$rn)
colnames(seq.abundancesV4_all)[colnames(seq.abundancesV4_all) == c("rn")] <- c("ASV")

taxa1.silva_V4 <- setDT(as.data.frame(taxa_V4_all), keep.rownames = TRUE)
colnames(taxa1.silva_V4)[colnames(taxa1.silva_V4) == c("rn")] <- c("ASV_seq")
taxa.abundance.silva_V4 <- merge(seq.abundancesV4_all, taxa1.silva_V4, by = "ASV_seq")
write.csv(taxa.abundance.silva_V4, "Output/dada2/taxa.silva_V4_all_boot95.csv", row.names = FALSE)

track_V4_all <- data.frame(cbind(rbind(track_V4_r1, track_V4_r2),c(rep("r1",nrow(track_V4_r1)), rep("r2",nrow(track_V4_r2)))))
track_V4_all <- cbind(rownames(track_V4_all), track_V4_all)
colnames(track_V4_all)[colnames(track_V4_all) == c("V5")] <- "run"

track_V4_all <- data.frame(sample = rownames(track_V4_all), track_V4_all)
write.csv(track_V4_all, file = "Output/dada2/track_reads_V4.csv", row.names = FALSE)

sessionInfo()