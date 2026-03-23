### Global ---------------------------------------------------------------------
set.seed(123)

source("Data_prep_models_vis.R")

library(glmmTMB)
library(visreg)
library(performance)
library(factoextra)
library(car)

### Diversity analysis between glaciers - models --------------------------
## V4
# Richness - Linear model
m_V4_R.1 <- glmmTMB(V4_observed ~ rad_total + rad_ratio:rad_total + Size_kmq + Average_altitude, 
                  family = "nbinom2", REML = TRUE, data = df_glac_dv); summary(m_V4_R.1)
check_model(m_V4_R.1); check_overdispersion(m_V4_R.1)

Anova(m_V4_R.1, type = "III", test = "Chisq")

# Shannon-W Index - Linear model
m_V4_S.1 <- glmmTMB(V4_Shannon ~ rad_total + rad_ratio:rad_total + Size_kmq + Average_altitude, 
               data = df_glac_dv, REML = TRUE); summary(m_V4_S.1)
check_model(m_V4_S.1)

Anova(m_V4_S.1, type = "III", test = "Chisq")

# Evenness - Linear model
m_V4_P.1 <- glmmTMB(V4_Pielou ~ rad_total + rad_ratio:rad_total + Size_kmq + Average_altitude, 
                    data = df_glac_dv, REML = TRUE); summary(m_V4_P.1)
check_model(m_V4_P.1)

Anova(m_V4_P.1, type = "III")

## V9
# Richness - Linear model
m_V9_R.1 <- glmmTMB(V9_observed ~ rad_total + rad_ratio:rad_total + Size_kmq + Average_altitude, 
                  REML = TRUE, family = "nbinom2", data = df_glac_dv); summary(m_V9_R.1)
check_model(m_V9_R.1); check_overdispersion(m_V9_R.1)

Anova(m_V9_R.1, type = "III")

### Animals dependence on radioactivity on glaciers -----------------------
df_glac_tar <- aggregate(. ~ glacier, df_sample_meta[,c(1, 4, 5)], FUN = sum)
df_glac_tar <- merge(df_glac_tar, aggregate(. ~ glacier, df_sample_meta[,c(1, 6, 8, 9)], FUN = mean), by = "glacier")

df_V9_tar <- subset(df_V9_s, Order == "Tardigrada")
df_V9_tar_num <- df_V9_tar[,3:159]

df_V9_tar_num <- df_V9_tar_num[,colSums(df_V9_tar_num) > 0]
df_V9_tar_ASV <- nrow(df_V9_tar_num) - colSums(df_V9_tar_num==0)

mean(df_V9_tar_ASV); sd(df_V9_tar_ASV)
table(df_V9_tar_ASV)

df_glac_tar$rad_total <- df_glac_tar$Cs137 + df_glac_tar$Pb210
df_glac_tar$rad_ratio <- df_glac_tar$Pb210/df_glac_tar$Cs137

mT.1 <- glmmTMB(tard_count ~ rad_total + rad_total:rad_ratio + offset(log(volume_ml)), 
                family = "nbinom2", REML = TRUE, data = df_glac_tar); summary(mT.1)

check_model(mT.1); check_overdispersion(mT.1)
visreg(mT.1, scale = "response")
Anova(mT.1, type = "III")

sessionInfo()