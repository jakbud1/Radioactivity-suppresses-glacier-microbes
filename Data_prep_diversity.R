### Global ----------------------------------------------------------------
set.seed(123)

source("ASV_post_processing.R")

library(phyloseq)
library(vegan)
library(tidyr)
library(readxl)

df_sample_meta <- read_xlsx("Input/meta_samples.xlsx") %>% 
  column_to_rownames("ID") 

### V4 - Diversity ----------------------------------------------
## Phyloseq objects
ASV_V4 <- otu_table(df_V4_asv, taxa_are_rows = TRUE)
TAX_V4 <- tax_table(df_V4_tax)
samples <- sample_data(df_sample_meta)

V4_ps <- phyloseq(ASV_V4, TAX_V4, samples)

sample_names(V4_ps); rank_names(V4_ps)

all(sample_names(V4_ps) %in% rownames(df_sample_meta))

## Rarefraction curve
tiff("Output/V4_rarecurve.tiff", units = "px", width = 2000, height = 1100, pointsize = 16)
rarecurve(t(df_V4_asv), step = 200, cex = 0.75, cex.axis = 2.2, cex.lab = 1.5 ,col = "#2C3930", ylab = "16S V4 rDNA ASVs", labelsize = 10)
dev.off()

## Normalisation 
V4_ps <- rarefy_even_depth(V4_ps, sample.size = 10000, rngseed = 1, replace = FALSE)

## Export diversity indexes 
V4_out_rich <- estimate_richness(V4_ps, measure = c("Observed", "Shannon"))

## Export evenness index
V4_evennes <- data.frame(ID = substring(rownames(V4_out_rich),2), Pielou = V4_out_rich$Shannon / log(V4_out_rich$Observed))

### V9 - Diversity ----------------------------------------------
## Phyloseq objects 
ASV_V9 <- otu_table(df_V9_asv, taxa_are_rows = TRUE)
TAX_V9 <- tax_table(df_V9_tax)
samples <- sample_data(df_sample_meta)

V9_ps <- phyloseq(ASV_V9, TAX_V9, samples); rm(ASV_V9, TAX_V9, df_V9_tax)

sample_names(V9_ps); rank_names(V9_ps)

all(sample_names(V9_ps) %in% rownames(df_sample_meta))
## Rarefraction curve
tiff("Output/V9_rarecurve.tiff", units = "px", width = 2000, height = 1100, pointsize = 16)
rarecurve(t(df_V9_asv),  step = 200, cex = 0.75, cex.axis = 2.2, cex.lab = 1.5, col = "#A27B5C", ylab = "18S V9 rDNA ASVs")
dev.off()

## Normalisation 
V9_ps <- rarefy_even_depth(V9_ps, rngseed = 1, replace = FALSE)

## Export diversity indexes
V9_out_rich <- estimate_richness(V9_ps, measure = c("Observed"))

### Calculate relative abundance averaged for each Phylum (V4) and Class (V9) -----------------------
## V4
V4_rlt_in <- df_V4_s

V4_samples_keep <- c(1,2,which(colnames(V4_rlt_in[,3:159]) %in% colnames(V4_ps@otu_table)) + 2,160:168)
V4_rlt_in <- V4_rlt_in[,V4_samples_keep]

V4_rlt_in$Phylum <- as.character(ifelse(is.na(df_V4_s$Phylum), "Unclassified", df_V4_s$Phylum))

Phyla_V4 <- unique(V4_rlt_in$Phylum)

V4_rlt <- data.frame(matrix(NA, nrow = 151, ncol = length(Phyla_V4) + 1))

for (i in 1:length(Phyla_V4)) {
  V4_rlt[,i] <- (colSums(subset(V4_rlt_in, Phylum == Phyla_V4[i])[,c(3:153)])/colSums(V4_rlt_in[,c(3:153)]))*100
  colnames(V4_rlt)[i] <- Phyla_V4[i]
}

V4_rlt[max(length(V4_rlt))] <- colnames(V4_rlt_in[c(3:153)])
colnames(V4_rlt)[max(length(V4_rlt))] <- "ID"

# Test whether this loop work properly:
table(rowSums(V4_rlt[,c(1:37)]))

## V9
V9_rlt_in <- df_V9_s
V9_samples_keep <- c(1,2,which(colnames(V9_rlt_in[,3:159]) %in% colnames(V9_ps@otu_table)) + 2,160:170)
V9_rlt_in <- V9_rlt_in[,V9_samples_keep]

V9_rlt_in$Class <- as.character(ifelse(is.na(df_V9_s$Class), "Unclassified", df_V9_s$Class))

Class_V9 <- unique(V9_rlt_in$Class)

V9_rlt <- data.frame(matrix(NA, nrow = 157, ncol = length(Class_V9) + 1))

for (i in 1:length(Class_V9)) {
  V9_rlt[,i] <- (colSums(subset(V9_rlt_in, Class == Class_V9[i])[,c(3:159)])/colSums(V9_rlt_in[,c(3:159)]))*100
  colnames(V9_rlt)[i] <- Class_V9[i] 
}

V9_rlt[max(length(V9_rlt))] <- colnames(V9_rlt_in[c(3:159)])
colnames(V9_rlt)[max(length(V9_rlt))] <- "ID"

# Tests whether this loop work properly:
table((colSums(subset(df_V9_s, Class == "Metazoa")[,c(3:159)])/colSums(df_V9_s[,c(3:159)]))*100 == V9_rlt$Metazoa)
table(rowSums(V9_rlt[,c(1:26)]))

### Clear -----------------------------------------------------------------
rm(V4_ps, V9_ps, samples, 
   i, Class_V9, Phyla_V4, V9_rlt_in, V4_rlt_in, V4_samples_keep, V9_samples_keep)

sessionInfo()