### Global ---------------------------------------------------------------------
set.seed(123)

library(tidyverse)
library(readxl)
library(tibble)

df_V4 <- read.csv("Output/dada2/taxa.silva_V4_all_boot95.csv")
df_V9 <- read.csv("Output/dada2/taxa.silva_V9_all_boot60.csv")

### Check reads counts and distribution ----------------------------------------
## V4 
hist(nchar(df_V4$ASV_seq))
data.frame(colSums(df_V4[, c(3:184)]))

## V9
hist(nchar(df_V9$ASV_seq))
data.frame(colSums(df_V9[, c(3:182)]))

### Merge runs -----------------------------------------------------------------
## V4
df_V4_s <- as.data.frame(t(df_V4[,-c(1,2, 185:190)]) %>%
  data.frame() %>%
  group_by(., id = gsub('\\..*', '', rownames(.))) %>%
  summarise_all(sum) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t())

df_V4_s <- cbind(df_V4[,c(1,2)], df_V4_s, df_V4[,185:190])

## V9
df_V9_s <- as.data.frame(t(df_V9[,-c(1,2, 183:191)]) %>%
data.frame() %>%
  group_by(., id = gsub('\\..*', '', rownames(.))) %>%
  summarise_all(sum) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t())

df_V9_s <- cbind(df_V9[,c(1,2)], df_V9_s, df_V9[,183:191])

## Rename columns
names(df_V4_s)[c(9:166)] <- substring(names(df_V4_s), first = 2, last = 4)[c(9:166)]
names(df_V9_s)[c(11:167)] <- substring(names(df_V9_s), first = 2, last = 4)[c(11:167)]

### Remove too short reads and with too low count ------------------------------
## V4
# Size selection
hist(nchar(df_V4_s$ASV_seq))
df_V4_s <- df_V4_s[nchar(df_V4_s$ASV_seq) > 200 & nchar(df_V4_s$ASV_seq) < 220,] # Keep ASVs in range of 200-220bp

# Remove samples with low reads
data.frame(colSums(df_V4_s[, c(3:166)]))
df_V4_s <- df_V4_s[,-c(127)] # Sample 129 removed (only 13 reads)

# Check ASV with reads = 0 (test)
table(rowSums(df_V4_s[,c(3:165)]) == 0)
df_V4_s <- df_V4_s[!rowSums(df_V4_s[,c(3:165)]) == 0,]

# V9
hist(nchar(df_V9_s$ASV_seq))
V9_rm <- as.data.frame(colSums(df_V9_s[, c(3:167)]))[,1] < 300
df_V9_s <- df_V9_s[,!V9_rm]; rm(V9_rm) # Remove samples with reads lower than 300. 

# Check ASV with reads = 0 (test)
table(rowSums(df_V9_s[,c(3:166)]) == 0)
df_V9_s <- df_V9_s[!rowSums(df_V9_s[,c(3:166)]) == 0,]
table(rowSums(df_V9_s[,c(3:166)]) == 0)

### Sum and remove ASV found in n.controls ---------------------------------------
df_track <- data.frame(Region = c("V4", "V9"), Before = c(NA, NA), After = c(NA, NA), Controls_reads = c(NA, NA))

## V4
df_V4_s$n.con.sum <- Reduce("+", df_V4_s[,c(3:8)]); df_track$Controls_reads[1] <- sum(df_V4_s$n.con.sum)
df_V4_s <- df_V4_s[,-c(3:8)]

df_V4_s$sum.ASV <- Reduce("+", df_V4_s[,c(3:159,166)]); df_track$Before[1] <- sum(df_V4_s$sum.ASV)
df_V4_s$n.con.share <- (df_V4_s$n.con.sum/df_V4_s$sum.ASV)*100
hist(df_V4_s$n.con.share, breaks = 200)

df_V4_s <- df_V4_s[df_V4_s$n.con.share <= 2,]# remove all ASV where reads in n. controls were under 1% except one Cyanobacteria ASV (1.8%) 
df_V4_s <- df_V4_s[df_V4_s$n.con.share <= 1 | df_V4_s$Phylum == "Cyanobacteria",]
df_track$After[1] <- sum(Reduce("+", df_V4_s[,c(3:159,166)]))
df_V4_s <- df_V4_s[complete.cases(df_V4_s$ASV),]

table(rowSums(df_V4_s[,c(3:159)]) == 0)
table(colSums(df_V4_s[,c(3:159)]) == 0)

## V9
df_V9_s$n.con.sum <- Reduce("+", df_V9_s[,c(3:9)]); df_track$Controls_reads[2] <- sum(df_V9_s$n.con.sum)
df_V9_s <- df_V9_s[,-c(3:9)]

df_V9_s$sum.ASV <- Reduce("+", df_V9_s[,c(3:159,168)]); df_track$Before[2] <- sum(df_V9_s$sum.ASV)
df_V9_s$n.con.share <- (df_V9_s$n.con.sum/df_V9_s$sum.ASV)*100
hist(df_V9_s$n.con.share, breaks = 200)

df_V9_s <- df_V9_s[df_V9_s$n.con.share < 1,]# remove all ASV where reads in controls were under 1%
df_track$After[2] <- sum(Reduce("+", df_V9_s[,c(3:159)]))

table(rowSums(df_V9_s[,c(3:159)]) == 0)
df_V9_s <- df_V9_s[!rowSums(df_V9_s[,c(3:159)]) == 0,]
table(colSums(df_V9_s[,c(3:159)]) == 0)

rm(df_V4, df_V9)

### Remove Fungi and other reads from V9 output know as false possitives
df_V9_s <- subset(df_V9_s, Class != "Fungi" | is.na(Class))
df_V9_s <- subset(df_V9_s, Order != "Mammalia" | is.na(Order))
df_V9_s <- subset(df_V9_s, Order != "Nematoda" | is.na(Order))
df_V9_s <- subset(df_V9_s, Order != "Arthropoda" | is.na(Order))

### Prepare ASV data to PhyloSeq format --------------------------------------------
## V4
rownames(df_V4_s) <- NULL

df_V4_tax <- df_V4_s[,c(2,160:165)] %>%
  column_to_rownames("ASV") 
df_V4_tax <- as.matrix(df_V4_tax)

df_V4_asv <- df_V4_s[,c(2,3:159)] %>%
  column_to_rownames("ASV") 
df_V4_asv <- as.matrix(df_V4_asv)

## V9 
rownames(df_V9_s) <- NULL

df_V9_tax <- df_V9_s[,c(2,160:165)] %>%
  column_to_rownames("ASV") 
df_V9_tax <- as.matrix(df_V9_tax)

df_V9_asv <- df_V9_s[,c(2,3:159)] %>%
  column_to_rownames("ASV") 
df_V9_asv <- as.matrix(df_V9_asv)

dev.off(dev.list()["RStudioGD"])

### Save data ---------------------------------------------------------------
write.csv(subset(df_V4_s, n.con.share == 0)[,c(1,2,160:168)], "Output/Bacterial_taxa_clean.csv")
write.csv(subset(df_V4_s, n.con.share > 0)[,c(1,2,160:168)], "Output/Bacterial_taxa_inNControl.csv")

write.csv(subset(df_V9_s, n.con.share == 0)[,c(1,2,160:170)], "Output/Eukaryote_taxa_clean.csv")
write.csv(subset(df_V9_s, n.con.share > 0)[,c(1,2,160:170)], "Output/Eukaryote_taxa_inNControl.csv")

write.csv(df_track, "Output/Track_controls_drop.csv")

sessionInfo()