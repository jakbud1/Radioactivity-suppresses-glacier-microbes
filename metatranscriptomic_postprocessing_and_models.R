### Global ----------------------------
library(KEGGREST)
library(tidyr)
library(ggpicrust2)
library(dplyr)
library(glmmTMB)
library(visreg)
library(performance)
library(purrr)
library(broom)
library(data.table)
library(forcats)
library(RColorBrewer)
library(patchwork)
library(car)

## Samples
samples <- c("DOS_24_02", "DOS_24_03", "DOS_24_06", "DOS_24_08", "DOS_24_09", 
                 "DOS_24_10", "DOS_24_11", "DOS_24_12", "DOS_24_15", "DOS_24_16")

## Pathways
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

meta_trans <- readxl::read_xlsx("Input/DOS_24_meta.xlsx")

## DNA_data --------------------------------------------------------------------
for (i in samples){
  df <- read.table(paste("Input/humann3/DNA/", i, "_genefamilies_uns_ko_cpm.tsv", sep = ""), col.names = c("K", "CPM"))[-c(1:2),]
  df <- ko2kegg_abundance(data = df)
  df$ko <- rownames(df)
  
  tmpnames <- df
  
  colnames(tmpnames)[1] <- i
  
  if(i == "DOS_24_02"){
    out <- df
  }
  out <- full_join(out, tmpnames, by = "ko")
  rm(tmpnames, df)
}

DNA <- out[,-1]
df_DNA_repair <- data.frame(t(subset(DNA, ko %in% pathways_df$code)))
colnames(df_DNA_repair) <- df_DNA_repair[1,]
df_DNA_repair <- df_DNA_repair[-1,]

df_DNA_repair <- cbind(rownames(df_DNA_repair), df_DNA_repair)
colnames(df_DNA_repair)[1] <- "sample_ID"

rm(out, i)

### RNA_data -------------------------------------------------------------------
for (i in samples){
  df <- read.table(paste("Input/humann3/RNA/", i, "_genefamilies_uns_ko_cpm.tsv", sep = ""), col.names = c("K", "CPM"))[-c(1:2),]
  df <- ko2kegg_abundance(data = df)
  df$ko <- rownames(df)
  
  tmpnames <- df
  
  colnames(tmpnames)[1] <- i
  
  if(i == samples[1]){
    out <- df
  }
  out <- full_join(out, tmpnames, by = "ko")
  rm(tmpnames, df)
}

RNA <- out[,-1]
df_RNA_repair <- data.frame(t(subset(RNA, ko %in% pathways_df$code)))
colnames(df_RNA_repair) <- df_RNA_repair[1,]
df_RNA_repair <- df_RNA_repair[-1,]

df_RNA_repair <- cbind(rownames(df_RNA_repair), df_RNA_repair)
colnames(df_RNA_repair)[1] <- "sample_ID"

### Merging --------------------------------------------------------------------
df_repair <- left_join(df_RNA_repair, df_DNA_repair,  by = "sample_ID", suffix = c(".R", ".D"))

df_repair <- data.frame(apply(df_repair[,-1], 2, as.numeric))
df_repair$sample_ID <- df_RNA_repair$sample_ID

df_OUT <- left_join(df_repair, meta_trans, by = "sample_ID")

### Models ---------------------------------------------------------------------
df_OUT <- subset(df_OUT, sample_ID %in% samples)
df_OUT$rad_total_k <- df_OUT$rad_total/1000
  
# Models over loop
coef_df <- data.frame()
full_summaries <- list()

for (p in pathways_df$code) {
  
  response_var <- paste0(p, ".R")
  offset_var   <- paste0(p, ".D")
  
  fmla <- as.formula(paste(response_var, "~ rad_total_k + offset(log(", offset_var, "))"))
  
  m <- glmmTMB(fmla, family = lognormal, REML = TRUE, data = df_OUT)
  
  full_summaries[[p]] <- summary(m)
  
  s <- summary(m)
  coefs <- as.data.frame(s$coefficients$cond)  # conditional model coefficients
  coefs$term <- rownames(coefs)
  coefs$pathway <- p
  colnames(coefs)[1:4] <- c("estimate", "std.error", "z.value", "p.value")
  
  coef_df <- bind_rows(coef_df, coefs)
}

# BH correction per predictor
coef_df <- coef_df %>%
  group_by(term) %>%
  mutate(p_adjust_BH = round(p.adjust(p.value, method = "BH"), 3)) %>%
  ungroup()

table(subset(coef_df, term !="(Intercept)")$p_adjust_BH < 0.05)

### Visualization - linear relationships----------------------------------------
## Data 
df_metatrans_l <- df_OUT %>%
  tidyr::pivot_longer(
    cols = matches("\\.[RD]$"),
    names_to = c("feature", ".value"),
    names_pattern = "^(.*)\\.([RD])$"
  )

names(df_metatrans_l)[names(df_metatrans_l) == "R"] <- "RNA_abundance"
names(df_metatrans_l)[names(df_metatrans_l) == "D"] <- "DNA_abundance"

df_metatrans_l$RNA_t_DNA_ratio <- df_metatrans_l$RNA_abundance/df_metatrans_l$DNA_abundance

colnames(df_metatrans_l)[which(colnames(df_metatrans_l) == "feature")] <- "code"

df_metatrans_l <- left_join(df_metatrans_l, pathways_df, by = "code")

priority_mechs <- c("DNA-repairing", "Metabolism & stress response")
sig_pathways <- c("ko00030", "ko00450", "ko03440", "ko02020", "ko00010")

df_plot_mech <- df_metatrans_l %>%
  mutate(
    alpha_val = ifelse(code %in% sig_pathways, 1, 0.2),
    line_type = ifelse(code %in% sig_pathways, "solid", "dashed"),
    mechanism_order = factor(
      mechanism,
      levels = c(priority_mechs,
                 setdiff(unique(mechanism), priority_mechs))
    )
  )

ordered_pathways <- df_plot_mech %>%
  arrange(mechanism_order) %>%
  pull(pathways_short) %>%
  unique()

df_plot_mech <- df_plot_mech %>%
  mutate(
    pathways_short_ordered = factor(
      pathways_short,
      levels = ordered_pathways
    )
  )

## Plot
g1_tt <- ggplot(df_plot_mech,
                aes(x = rad_total,
                    y = RNA_t_DNA_ratio,
                    color = mechanism)) + 
  
  geom_point(aes(alpha = alpha_val), size = 2) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(group = code, linetype = line_type),
    size = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~pathways_short_ordered,
             scales = "free",
             nrow = 3) +
  scale_alpha_identity() +
  scale_linetype_manual(values = c("solid" = "solid",
                                   "dashed" = "dashed")) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  scale_y_log10(breaks = scales::log_breaks(n = 6)) +
  labs(
    x = expression("Total radioactivity (Bq kg"^-1*")"),
    y = "RNA transcripts per DNA copy",
    color = "Mechanism"
  ); g1_tt
# ggsave("Output/Fig. 4.png", g1_tt, dpi = 300, width = 12, height = 10, units = "in", bg = "white")

### Visualization - coefs ----------------------------------------
coef_df_m <- coef_df %>%
  subset(term == "rad_total_k") %>%
  mutate(
    signif = p_adjust_BH < 0.05,
    lower = estimate - std.error,
    upper = estimate + std.error
  )

colnames(coef_df_m)[6] <- "code"
coef_df_m <- left_join(coef_df_m, pathways_df, by = "code")

gg_coef <- ggplot(coef_df_m, 
       aes(x = estimate,
           y = reorder(pathways_short, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper,
                     color = mechanism,
                     alpha = signif),
                 height = 0.2,
                 size = 1) +
  geom_point(aes(color = mechanism,
                 size = signif,
                 alpha = signif)) +
  scale_alpha_manual(values = c(`TRUE` = 1,
                                `FALSE` = 0.3),
                     guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3,
                               `FALSE` = 2),
                    guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") + 
  labs(
    x = "Estimate (± SE)",
    y = "KEGG pathway"); gg_coef

### Combined plot ------------------------------------
combined_plot_meta <- (g1_tt | gg_coef) +
  plot_layout(widths = c(3, 1)) + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "", 
                  tag_suffix = ")", 
                  theme = theme(plot.tag = element_text(size = 14, face = "bold"))); combined_plot_meta

ggsave("Output/Fig. 5.png", combined_plot_meta, dpi = 300, width = 16, height = 8, units = "in", bg = "white")

### Models check for significant pathways -----------------
## Glycolysis
m.Gl <- glmmTMB(ko00010.R ~ rad_total_k + offset(log(ko00010.D)),family = lognormal, 
                                           REML = TRUE, data = df_OUT); summary(m.Gl)
check_model(m.Gl)
Anova(m.Gl)
           
## Se metabolism 
m.Se <- glmmTMB(ko00450.R ~ rad_total_k + offset(log(ko00450.D)),family = lognormal, 
                REML = TRUE, data = df_OUT); summary(m.Se)
check_model(m.Se)
Anova(m.Se)

## Two-component 
m.Tc <- glmmTMB(ko02020.R ~ rad_total_k + offset(log(ko02020.D)),family = lognormal, 
                REML = TRUE, data = df_OUT); summary(m.Tc)
check_model(m.Tc)
Anova(m.Tc)

## PPP
m.ppp <- glmmTMB(ko00030.R ~ rad_total_k + offset(log(ko00030.D)),family = lognormal, 
                REML = TRUE, data = df_OUT); summary(m.ppp)
check_model(m.ppp)
Anova(m.ppp)

## HR 
m.hr <- glmmTMB(ko03440.R ~ rad_total_k + offset(log(ko03440.D)),family = lognormal, 
                 REML = TRUE, data = df_OUT); summary(m.hr)
check_model(m.hr)
Anova(m.hr)
