### Global ----------------------------------------------------------------
source("Models.R")

library(RColorBrewer)
library(ggplot2)
library(data.table)
library(corrplot)
library(visreg)

col_V4 <- c("#D9D9D9", "#66C2A5", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4",  "#E6F5C9",
  "#FFF2AE", "#F1E2CC", "#CCEBC5", "#BC80BD", "#FFFFFF")
col_V9 <- c("#D9D9D9", "#66C2A5", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",  
  "#FFF2AE", "#FFFFFF")

### Add Tardigrada data to ASVs data --------------------------------------
df_animals_gl <- aggregate(. ~ glacier, df_sample_meta[,c(1, 4, 5)], FUN = sum)
df_animals_gl <- merge(df_animals_gl, df_glac_dv[,c(1:8,11,12)], by = "glacier")

df_animals_gl$Tard_dens <- df_animals_gl$tard_count/df_animals_gl$volume_ml

V4_cor <- merge(V4_rlt_wide[,c(1,3:12, 14:17)], df_animals_gl[,c(1,13)], by = "glacier" )
V9_cor <- merge(V9_rlt_wide[,c(1,3:8, 11:14)], df_animals_gl[,c(1,13)], by = "glacier" )

### Descriptive statistics ------------------------------------------------
## Radioactivity - for samples which were included in the diversity analysis, for all collected samples the results are in df_glac_tar df and they are reported in MS - results.
min(df_glac_dv$Pb210); max(df_glac_dv$Pb210)
mean(df_glac_dv$Pb210); sd(df_glac_dv$Pb210)
median(df_glac_dv$Pb210); IQR(df_glac_dv$Pb210)

min(df_glac_dv$Cs137); max(df_glac_dv$Cs137)
mean(df_glac_dv$Cs137); sd(df_glac_dv$Cs137)
median(df_glac_dv$Cs137); IQR(df_glac_dv$Cs137)

sd(subset(df_sample_meta, glacier == "Ventina")$Cs137)
sd(subset(df_sample_meta, glacier == "Mont Miné")$Pb210)

## Diversity
cor(df_animals_gl$V4_observed, df_animals_gl$V9_observed)

nrow(df_V4_s); ncol(df_V4_s) - 11 # Number V4 ASVs and samples at the end 
nrow(df_V9_s); ncol(df_V9_s) - 13 # Number V9 ASVs and samples at the end 

table(nchar(df_V4_s$ASV_seq)); median(nchar(df_V4_s$ASV_seq)) # V4 ASV seq length distribution
table(nchar(df_V9_s$ASV_seq)); median(nchar(df_V9_s$ASV_seq)) # V9 ASV seq length distribution

## Tardigrada
hist(df_glac_tar$tard_dens)
min(df_glac_tar$tard_dens); max(df_glac_tar$tard_dens)
median(df_glac_tar$tard_dens); IQR(df_glac_tar$tard_dens)

### Radioactivity data visualisation --------------------------------------
df_total <- df_glac_tar[,c(1,5,6)]
df_total$total_rad <- df_total$Pb210 + df_total$Cs137

df_box <- melt(as.data.table(df_sample_meta[,c(1,3,8,9)]), id.vars = c("code", "glacier"), variable.name = "Isotope")
df_box$glacier <- factor(df_box$glacier, levels = df_total[order(df_total$total_rad), ]$glacier)

levels(df_box$Isotope) <- c("Pb-210", "Cs-137")

g1 <- ggplot(df_box, aes(glacier, value, color = Isotope, fill = Isotope))

g1.o <- g1 + geom_boxplot(alpha = 0.5) +
  theme_classic(base_size = 14) + ylab(expression(Activity~concentration~(Bq~kg^{-1}))) + 
  scale_color_manual(values = c("#8FBC8F", "#E69F00")) + scale_fill_manual(values = c("#8FBC8F", "#E69F00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(),  legend.position = "right",
        legend.justification = "left",
        legend.box.just = "left") + xlab("Glacier");g1.o
# ggsave("Output/Fig. 1.tiff", dpi = 300, scale = 1, width = 2900, height = 1500, units = "px", bg = "white")

length(unique(df_box[complete.cases(df_box$value),]$code)) # how many samples used in analysis?

### Diversity - Relative share variation between glaciers -----------------
# V4 - relative share 
V4_rlt_long_gg <- V4_rlt_long
V4_rlt_long_gg$glacier <- factor(V4_rlt_long_gg$glacier, levels = df_glac_dv[order(df_glac_dv$rad_total), ]$glacier)

gg_div_g_V4 <- ggplot(V4_rlt_long_gg, aes(x = glacier, y = Relative_share, fill = Phylum)) + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 14) + 
  ylab("Relative abundance of bacteria (%)") + scale_fill_manual(values = col_V4) + 
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.text = element_text(size = 14, colour = "black"),  legend.position = "right",
        legend.justification = "left",
        legend.box.just = "left", legend.key.size = unit(0.7,"line"), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)); gg_div_g_V4
# ggsave("Output/Fig. 2a.tiff", dpi = 300, scale = 1, width = 3000, height = 1600, units = "px", bg = "white")

# V9 - relative share
V9_rlt_long_gg <- V9_rlt_long
V9_rlt_long_gg$glacier <- factor(V9_rlt_long_gg$glacier, levels = df_glac_dv[order(df_glac_dv$rad_total), ]$glacier)

gg_div_g_V9 <- ggplot(V9_rlt_long_gg, aes(x = glacier, y = Relative_share, fill = Class)) + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 14) + 
  ylab("Relative abundance of eukaryotes (%)") + xlab("Glacier") + scale_fill_manual(values = col_V9) + 
  guides(fill=guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.text = element_text(size = 14, colour = "black"),  legend.position = "right",
        legend.justification = "center",
        legend.box.just = "center", legend.margin = margin(), legend.key.size = unit(0.7,"line"), axis.text.x = element_text(angle = 45, hjust = 1)); gg_div_g_V9
# ggsave("Output/Fig. 2b.tiff", dpi = 300, scale = 1, width = 3000, height = 1600, units = "px", bg = "white")

### Diversity - correlations ----------------------------------------------
## V4 
names(V4_cor)[12:16] <- c("Pb-210", "Cs-137", 
                          "Glacier altitude", "Glacier size", "Tardigrada density")

V4_cor_m <- cor(V4_cor[,-1])

png("Output/S1.a.png", width = 1500, height = 1600, pointsize = 40)
corrplot(V4_cor_m, method = 'circle', type = 'lower', insig = 'blank',
         number.cex = 0.8, diag = FALSE, tl.col = 'black', tl.srt = 45)
dev.off()

## V9 
names(V9_cor)[8:12] <- c("Pb-210", "Cs-137", 
                          "Glacier altitude", "Glacier size", "Tardigrada density")

V9_cor_m <- cor(V9_cor[,-1])

png("Output/S1.b.png", width = 1500, height = 1600, pointsize = 40)
corrplot(V9_cor_m, method = 'circle', type = 'lower', insig = 'blank',
         number.cex = 0.8, diag = FALSE, tl.col = 'black', tl.srt = 45)
dev.off()

### Models results --------------------------------------------------------
## V4 
vis.V4.1 <- visreg(m_V4_R.1, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)
vis.V4.1.g <- ggplot(vis.V4.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#9E9D24") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.15, col = "white") + 
  geom_point(data = data.frame(m_V4_R.1$frame), aes(x = rad_total, y = V4_observed), shape = 2, col = "#9E9D24") +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Average number of ASV on glacier"); vis.V4.1.g 
ggsave("Output/Fig. 2a.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

vis.V4.2 <- visreg2d(m_V4_R.1, "rad_total", "rad_ratio", scale = "response", plot.type = "gg", 
                     col = c(
                       "#543005", "#7F511E",  "#A6611A",  "#D07C32",  "#C2A83E",  "#9E9D24", 
                       "#BFCB78", "#E6E8A0", "#FCFCD6")) + theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab(expression(paste(""^"210", "Pb", "/"^"137", "Cs (activity ratio)"))) + 
  theme(legend.position = "bottom")
ggsave("Output/Fig. 2c.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

## V9 
vis.V9.1 <- visreg(m_V9_R.1, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)
ggplot(vis.V9.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#009999") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.15, col = "white") + 
  geom_point(data = data.frame(df_glac_dv), 
             aes(x = rad_total, y = V9_observed), color = "#009999", shape = 21) +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Average number of ASV on glacier")
ggsave("Output/Fig. 2b.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

vis.V9.2 <- visreg2d(m_V9_R.1, "rad_total", "rad_ratio", scale = "response", plot.type = "gg", 
                     col = c(
                       "#543005", "#7F511E", "#A6611A", "#C08C46", "#7FBF9D", "#40B5AD", "#009999", "#7FD3E0",  
                       "#E0F7FA"))
vis.V9.2 + theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab(expression(paste(""^"210", "Pb", "/"^"137", "Cs (activity ratio)"))) + 
  theme(legend.position = "bottom")
ggsave("Output/Fig. 2d.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

## Tardigrada 
vis.tar.1 <- visreg(mT.1, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)

ggplot(vis.tar.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#d65917") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.15, col = "white") + 
  geom_point(data = data.frame(mT.1$frame), aes(x = rad_total, y = tard_count), shape = 21, col = "#d65917") +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Tardigrada abundance")
ggsave("Output/Fig. 2e.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")


### Combined plots ---------------------- 
## Radioactivity and community composition
combined_plot_1 <- g1.o / gg_div_g_V4 / gg_div_g_V9; combined_plot_1
ggsave("Output/Fig. 1.png", combined_plot_1, dpi = 300, width = 10, height = 12, units = "in", bg = "white", scale = 0.9)
