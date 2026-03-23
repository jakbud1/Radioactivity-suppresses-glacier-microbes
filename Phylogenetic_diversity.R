### Global ---------------------------------------------------------------------
set.seed(123)

library(tidyverse)
library(phyloseq)
library(ape)
library(DECIPHER)
library(phangorn)
library(vegan)
library(picante)
library(glmmTMB)
library(performance)
library(car)
library(ggplot2)

### Data ----------------------------------------------------------------------
asv_table_PD <- subset(df_V4_s, Kingdom == "Bacteria")[,-c(160:168)]

# ASV seq 
asv_seqs_PD <- asv_table_PD[,c(1:2)]

# ASV abundance
rownames(asv_table_PD) <- asv_table_PD$ASV
asv_table_PD <- asv_table_PD[,-c(1:2)]

# metadata 
asv_samples_PD <- df_cryoconite_dv[,c(1:2,9,10)]

# Checks
all(asv_table_PD %% 1 == 0)
any(asv_table_PD < 0)
table(colSums(asv_table_PD) == 0)

# Rare ASV  removing
asv_table_PD <- asv_table_PD[rowSums(asv_table_PD) >= 10, ]

# Sequence table synq.
asv_seqs_PD <- asv_seqs_PD[asv_seqs_PD$ASV %in% rownames(asv_table_PD), ]

### Phyloseq prep --------------------------------------------------------------
# OTU table
otu_phy <- otu_table(as.matrix(asv_table_PD), taxa_are_rows = TRUE)

# Sample metadata
metadata_phy <- asv_samples_PD %>%
  column_to_rownames("ID") %>%
  select(glacier, Pb210, Cs137)
sample_phy <- sample_data(metadata_phy)

### Phylogenetic tree from ASV -------------------------------------------------
dna <- DNAStringSet(asv_seqs_PD$ASV_seq)
names(dna) <- asv_seqs_PD$ASV

# Aligning sequences
length(dna) ### ---------- check
nrow(asv_table_PD)

alignment <- AlignSeqs(dna, anchor=NA)

alignment_matrix <- as.matrix(alignment)

# Alignment to phyDat
phydat <- phyDat(alignment_matrix, type = "DNA")

# Distance matrix
dm <- dist.ml(phydat)

# Build tree
treeNJ <- NJ(dm)
treeNJ <- ladderize(treeNJ)

### Rarefaction ----------------------------------------------------------------
sample_depths <- colSums(asv_table_PD)
summary(sample_depths)

asv_rare_PD <- t(rrarefy(t(asv_table_PD), sample = 10000))

### Match tree and ASV data ----------------------------------------------------
asvs_in_data <- rownames(asv_rare_PD)
phy_tree <- drop.tip(
  treeNJ,
  setdiff(treeNJ$tip.label, rownames(asv_rare_PD))
)

### Faith’s PD -----------------------------------------------------------------
faith_pd <- pd(t(asv_rare_PD), phy_tree, include.root = FALSE)

metadata_phy$sample_id <- as.character(rownames(metadata_phy))
pd_df <- data.frame(
  sample_id = rownames(faith_pd),
  PD = faith_pd$PD,
  SR = faith_pd$SR
) %>%
  left_join(metadata_phy, by = "sample_id")

pd_df <- pd_df[complete.cases(pd_df$Pb210),]

pd_glacier <- pd_df %>%
  group_by(glacier) %>%
  summarise(
    PD_mean = mean(PD),
    PD_sd   = sd(PD),
    SR_mean = mean(SR),
    Pb210_mean = mean(Pb210),
    Cs137_mean = mean(Cs137),
    n = n()
  )

pd_glacier$rad_total <- pd_glacier$Pb210_mean + pd_glacier$Cs137_mean
pd_glacier$rad_ratio <- pd_glacier$Pb210_mean/pd_glacier$Cs137_mean

### Models ---------------------------------------------------------------------
m_pd <- glmmTMB(PD_mean ~ rad_total + rad_total:rad_ratio, 
                data = pd_glacier, REML = TRUE); summary(m_pd)
check_model(m_pd)
Anova(m_pd)

### Visualisation --------------------------------------------------------------
PD_g <- ggplot(pd_glacier, aes(x = rad_total, y = PD_mean)) +
  geom_point(shape = 2, color = "#C2A83E") +
  geom_smooth(method = "lm", se = TRUE, color = "#C2A83E", alpha = 0.3, linewidth = 0.5) +
  ylab("Faith’s Phylogenetic Diversity") + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) +
  theme_classic(); PD_g

ggsave("Output/Fig. 2f.png", PD_g, scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

sessionInfo()