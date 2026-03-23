### For PICRUSt2, V4 16S rDNA (Bacteria) data were used.
library(biomformat)
library(seqinr)
library(spgs)

### Preprocess ASV data
## data 
df_V4_s_nc <- df_V4_s
df_V4_s_nc$Phylum <- ifelse(is.na(df_V4_s_nc$Phylum), "NA", df_V4_s_nc$Phylum)
df_V4_s_nc <- subset(df_V4_s_nc, Phylum != "Cyanobacteria") # Cyanobacteria are low-diversity, yet they can possess more abundant pathways to withstand high UV stress, and are also known to be positively correlated with radioisotope content. To reduce false positives, we removed them from the PICRUSt2 inference.

## Remove extreme low-abundant taxa 
df_V4_s_nc_cut <- df_V4_s_nc[rowSums(df_V4_s_nc[,3:159]) > 5,]

## .biom file
df_abu_nc <- df_V4_s_nc_cut[,2:159]
row.names(df_abu_nc) <- df_abu_nc$ASV
df_abu_nc <- df_abu_nc[,-1]

ASV_biom_nc <- make_biom(data = df_abu_nc)
write_biom(ASV_biom_nc, "Output/dada2/V4_ASV_abun_nc.biom")

## .fna file
write.fasta(as.list(reverseComplement(df_V4_s_nc_cut$ASV_seq)), as.list(df_V4_s_nc_cut$ASV), 
            "Output/dada2/V4_ASV_seqs_nc.fna")
