### Global ---------------------------------------------------------------------
set.seed(123)

library(ggpicrust2)
library(tidyverse)
library(patchwork)
library(stringr)
library(vegan)

### Input ----------------------------------------------------------------------
ko_abundance_t_nc <- read_delim("Input/Picrust2/pred_metagenome_unstrat_nc.tsv", delim = "\t")

metadata_t <- tibble::tibble(df_cryoconite_dv[,1:2])
colnames(metadata_t)[1] <- "sample_name"

### Data preparation -----------------------------------------------------------
glaciers_class <- df_glac_dv[,c(1, 13)]
glaciers_class <- glaciers_class[order(glaciers_class$rad_total),]
glaciers_class <- subset(glaciers_class, glacier != "Mont Miné")

glaciers_class$Contamination <- c(rep("Low", 7), rep("Mid", 1), rep("High", 7))

metadata <- subset(metadata_t, glacier %in% 
                        subset(glaciers_class, Contamination != "Mid")$glacier) 
metadata <- left_join(metadata, glaciers_class, by = "glacier")[,-3]

metadata_avg <- metadata %>%
  select(glacier, Contamination) %>%   
  distinct()                            

ko_abundance_nc <- ko_abundance_t_nc[, c("function", metadata$sample_name)]

colnames(ko_abundance_nc)[-1] <- metadata$glacier
colnames(ko_abundance_nc)[1] <- "temp"

ko_abundance_avg <- cbind(
  temp = ko_abundance_nc$temp,
  as.data.frame(
    sapply(
      split.default(
        ko_abundance_nc[ , !names(ko_abundance_nc) %in% "temp"],
        names(ko_abundance_nc)[!names(ko_abundance_nc) %in% "temp"]
      ),
      rowMeans,
      na.rm = TRUE
    )
  )
)

colnames(ko_abundance_avg)[1] <- "function"

# asses mean values for groups
mean(subset(metadata_phy, glacier %in% metadata_avg$glacier)$Pb210, na.rm = TRUE) ### The metadata_phy is first defined in Phylogenetic_diversity.R

mean_ac <- full_join(metadata_phy, metadata_avg, by = "glacier")

mean(subset(mean_ac, Contamination == "Low")$Pb210, na.rm = TRUE)
mean(subset(mean_ac, Contamination == "Low")$Cs137, na.rm = TRUE)

mean(subset(mean_ac, Contamination == "High")$Pb210, na.rm = TRUE)
mean(subset(mean_ac, Contamination == "High")$Cs137, na.rm = TRUE)

# ko to kegg (pathways) abundance
kegg_abundance_avg <- ko2kegg_abundance(data = ko_abundance_avg)

kegg_abundance_avg <- kegg_abundance_avg[complete.cases(kegg_abundance_avg[1]),]

rm(metadata_t)

# Selecting pathways related to DNA-repair and oxidative stress systems 
pathways_df <- tibble::tibble(
  pathways = c(
    "Base excision repair",
    "Nucleotide excision repair",
    "Mismatch repair",
    "Homologous recombination",
    "Non-homologous end joining",
    "Pentose phosphate pathway (NADPH)",
    "Glutathione metabolism",
    "Selenocompound metabolism",
    "Protein folding",
    "Protein processing",
    "ABC transporters",
    "Two-component system",
    "DNA replication",
    "Bacterial translation",
    "Central carbon metabolism"
  ),
  code = c(
    "ko03410",
    "ko03420",
    "ko03430",
    "ko03440",
    "ko03450",
    "ko00030",
    "ko00480",
    "ko00450",
    "ko04141",
    "ko04142",
    "ko02010",
    "ko02020",
    "ko03030",
    "ko03010",
    "ko00010"
  ), 
  mechanism = c(
    "DNA repair",
    "DNA repair",
    "DNA repair",
    "DNA repair",
    "DNA repair",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response",
    "Metabolism & stress response"
  ), 
  pathways_short = c(
    "BER",
    "NER",
    "MMR",
    "HR",
    "NHEJ",
    "PPP",
    "GSH",
    "Se metabolism",
    "Protein folding",
    "Protein processing",
    "ABC transporters",
    "Two-component",
    "DNA replication",
    "Translation",
    "Glycolysis"
  )
)

kegg_abundance_avg_sel <- kegg_abundance_avg[rownames(kegg_abundance_avg) %in% pathways_df$code, ]

### Differential abundance analysis based on ALDEx2 ----------------------------
daa_results_avg <- pathway_daa(abundance = kegg_abundance_avg_sel, metadata = metadata_avg,
                                    group = "Contamination", daa_method = "ALDEx2", select = NULL, p.adjust = "BH", 
                                   reference = "Low", include_abundance_stats = TRUE)

daa_results_avg$p_adjust <- round(daa_results_avg$p_adjust, 3)

daa_results_avg <- subset(daa_results_avg, method == "ALDEx2_Welch's t test")

### Visualization of pathways change in relative abundance ---------------------
# Convert to long format
df_long <- daa_results_avg %>%
  pivot_longer(
    cols = starts_with("mean_rel_abundance_group"),
    names_to = "Group",
    values_to = "Mean"
  ) %>%
  mutate(
    SD = case_when(
      Group == "mean_rel_abundance_group1" ~ sd_rel_abundance_group1,
      Group == "mean_rel_abundance_group2" ~ sd_rel_abundance_group2
    ),
    Group = case_when(
      Group == "mean_rel_abundance_group1" ~ "High",
      Group == "mean_rel_abundance_group2" ~ "Low"
    ),
    Significant = ifelse(p_adjust < 0.05, "Yes", "No")
  )

# Merge descriptive names
df_long_named <- df_long %>%
  left_join(pathways_df, by = c("feature" = "code"))

# Plot
g1_pathways <- ggplot(df_long_named, aes(x = reorder(pathways_short, -log2_fold_change), y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.3,
                position = position_dodge(width = 0.8)) +
  facet_wrap(.~mechanism, scales = "free") +
  scale_fill_manual(values = c("High" = "#e31a1c", "Low" = "#1f78b4")) +
  geom_text(data = df_long_named %>% filter(Group == "High"),
            aes(x = pathways_short, y = Mean + 0.001, label = ifelse(Significant == "Yes", "*", "")),
            position = position_dodge(width = 0.8), vjust = 0, size = 8) +
  labs(
    x = "KEGG pathway",
    y = "Mean relative abundance",
    fill = "Radioisotope activity concentrations"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top"); g1_pathways

### Visualization of pathways change in fold-change ---------------------
# Data prepare for fold-change barplot
daa_results_avg_t <- daa_results_avg %>%
  mutate(Significant = ifelse(p_adjust < 0.05, "Yes", "No"))
daa_results_avg_t$log2_fold_change <- -daa_results_avg_t$log2_fold_change

colnames(daa_results_avg_t)[1] <- "code"
daa_results_avg_t <- left_join(daa_results_avg_t, pathways_df, by = "code")

# Fold-change barplot
fc_plot <- ggplot(daa_results_avg_t, 
                  aes(x = reorder(pathways_short, log2_fold_change), y = log2_fold_change, fill = Significant)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(.~mechanism, scales = "free_x") + 
  scale_fill_manual(values = c("Yes" = "black", "No" = "grey70")) +
  geom_text(aes(label = ifelse(Significant == "Yes", "*", "")),
            nudge_y = 0.02,  
            size = 8,     
            vjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +  # add space above bars
  labs(x = "KEGG pathway", y = "Log2 fold change (high / low)", fill = "Significant") +
  theme_classic(base_size = 14) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)); fc_plot

# Volcano plot
volcano_plot <- ggplot(daa_results_avg_t, aes(x = log2_fold_change, y = -log10(p_adjust), color = Significant)) +
  geom_point(size = 3) + geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_color_manual(values = c("Yes" = "black", "No" = "grey50")) +
  geom_text(data = daa_results_avg_t %>% filter(Significant == "Yes"),
            aes(label = "NHEJ"), vjust = 2, hjust = 0.7, size = 5) + guides(fill = "none", color = "none") + 
  labs(x = "Log2 fold change (high / low)", y = "-Log10 adjusted p-value", color = "Significant") + 
  theme(axis.title.x = element_text(vjust = 1, margin = margin(t = 3))) + 
  theme_classic(base_size = 14); volcano_plot

### Functional abundances variation between High and Low radioactive groups ----
## Data 
ko_pa <- ko_abundance_avg[,-1]
ko_pa[ko_pa > 0] <- 1

func_richness <- colSums(ko_pa)

rich_df <- data.frame(
  Glacier = names(func_richness),
  Func_Richness = as.numeric(func_richness)
)

rich_df <- merge(rich_df, metadata_avg, by.x = "Glacier", by.y = "glacier")

## Test: Fun. richness ~ radioactivity status
wilcox.test(Func_Richness ~ Contamination, data = rich_df)

## Test: Fun. Dissimilarity ~ radioactivity status
bray_pa <- vegdist(t(ko_abundance_avg[,-1]), method = "bray")
adonis2(bray_pa ~ Contamination, data = metadata_avg)

## Visualization of raw functional richness 
g_fun_d <- ggplot(rich_df, aes(Contamination, Func_Richness, fill = Contamination)) +  
geom_boxplot(alpha = 0.5) + geom_jitter(alpha = 0.5) + theme_classic(base_size = 14) + 
  scale_fill_manual(values = c("High" = "#e31a1c", "Low" = "#1f78b4")) + 
  ylab("Functionall richness \n (KEGG Orthologs)") + 
  xlab("Radioisotope activity concentrations") + theme(legend.position = "none");g_fun_d

## Visualization of changes in abundance of redundant functional orthologs 
# Data 
high_samples <- metadata_avg$glacier[metadata_avg$Contamination == "High"]
low_samples  <- metadata_avg$glacier[metadata_avg$Contamination == "Low"]

mean_low_all  <- rowMeans(ko_abundance_avg[, low_samples])

# logical vector: TRUE if lost in High
lost_flag <- (rowSums(ko_pa[, high_samples]) == 0 &
                rowSums(ko_pa[, low_samples]) > 0)

# values for the two groups
mean_low_lost      <- mean_low_all[lost_flag]
mean_low_retained  <- mean_low_all[!lost_flag]

plot_df <- data.frame(
  mean_low = c(mean_low_lost, mean_low_retained),
  status = c(rep("Lost", length(mean_low_lost)),
             rep("Retained", length(mean_low_retained)))
)

# Boxplots of lost functions abundance 
g_fun_a <- ggplot(plot_df, aes(x = status, y = mean_low)) +
  geom_jitter(alpha = 0.1, size = 0.7) + geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "",
    x = "Functional response on glaciers with high \n radiositopes activity concentration",
    y = "KEGG Orthologs abundance on glaciers with low \n radioisotopes activity concentration (log scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)); g_fun_a

### Combined plots -------------------------------------------------------------
## Functional richness and redundancy
combined_plot_fun <- (g_fun_d | g_fun_a) + 
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  ); combined_plot_fun
ggsave("Output/Fig. 3.png", combined_plot_fun, dpi = 300, width = 10, height = 5, units = "in", bg = "white")

## Pathways abundance and log-change 
combined_plot_path <- 
  g1_pathways /
  ((fc_plot | volcano_plot) + plot_layout(widths = c(0.65, 0.35))) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "", 
                  tag_suffix = ")", 
                  theme = theme(plot.tag = element_text(size = 14, face = "bold"))); combined_plot_path
ggsave("Output/Fig. 4.png", combined_plot_path, dpi = 300, width = 10, height = 10, units = "in", bg = "white")

### Save results ---------------------------------------------------------------
write.csv(kegg_abundance_avg_sel, "Output/Kegg_abundance_avg_sel.csv")
write.csv(daa_results_avg, "Output/DAA_results.csv")

sessionInfo()