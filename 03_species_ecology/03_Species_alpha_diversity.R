###############################################################################
### LOAD LIBRARIES
###############################################################################
library(tidyverse)
library(vegan)
library(FSA)           
library(ggpubr)
library(rstatix)
library(sandwich)      
library(lmtest)
library(broom)

###############################################################################
### 1. LOAD DATA
###############################################################################
species_mat <- read.csv("TPM_raw_filtered_bact.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_metadata_aligned_bact.csv")

###############################################################################
### 2. PREPARE SAMPLE × SPECIES MATRIX & REMOVE BAD METADATA
###############################################################################
species_mat_t <- as.data.frame(t(species_mat))
species_mat_t$Sample_ID <- rownames(species_mat_t)

# Remove samples missing metadata
meta <- meta %>% drop_na(Group, Age, BMI, Sex)

# Keep species TPM for only samples in metadata
A <- species_mat_t %>%
  filter(Sample_ID %in% meta$Sample_ID) %>%
  left_join(meta, by="Sample_ID")

###############################################################################
### 3. PREPARE NUMERIC ABUNDANCE MATRIX
###############################################################################
abund <- A %>% select(-Sample_ID, -Group, -BMI, -Age, -Sex)

abund <- abund %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))
abund[is.na(abund)] <- 0

###############################################################################
### 4. COMPUTE ALPHA DIVERSITY METRICS
###############################################################################
A$Shannon  <- diversity(abund, index = "shannon")
A$Simpson  <- diversity(abund, index = "simpson")
A$Richness <- rowSums(abund > 0)
A$Pielou   <- ifelse(A$Richness > 1, A$Shannon / log(A$Richness), NA)

cat("Alpha diversity metrics computed\n")

###############################################################################
### 5. PUBLICATION-QUALITY ALPHA DIVERSITY BOXPLOTS
###############################################################################

metrics <- c("Shannon", "Simpson", "Pielou", "Richness")
metric_labels <- c(
  Shannon = "Shannon Diversity",
  Simpson = "Simpson Diversity",
  Pielou  = "Pielou Evenness",
  Richness = "Species Richness"
)

group_colors <- c(
  Healthy="#ffcc00",
  Adenoma="#984ea3",
  Cancer ="#4daf4a"
)

for (m in metrics) {
  
  A$MetricValue <- A[[m]]
  
  p <- ggplot(A, aes(Group, MetricValue)) +
    geom_jitter(width = 0.25, alpha = 0.1, size = 2, color="black") +
    geom_boxplot(
      aes(fill = Group),
      color="black", linewidth=1.25,
      width=0.55, alpha=0.85, outlier.shape=NA
    ) +
    scale_fill_manual(values = group_colors) +
    labs(x="", y=metric_labels[m], title=metric_labels[m]) +
    theme_bw(base_size=16) +
    theme(
      plot.title = element_text(face="bold", size=18, hjust=0.5),
      axis.title.y = element_text(face="bold", size=16),
      axis.text.x = element_text(angle=45, hjust=1, face="bold", size=14),
      axis.text.y = element_text(size=14),
      panel.grid = element_blank(),
      panel.border = element_rect(color="black", fill=NA, linewidth=1.3),
      legend.position="none"
    )
  
  print(p)
  
  ggsave(
    paste0("Species_", m, "_Alpha_Boxplot.png"),
    p, width=4, height=5, dpi=600
  )
}

###############################################################################
### 6. UNIFIED KRUSKAL–WALLIS + DUNN POST-HOC
###############################################################################
KW_all   <- tibble()
DUNN_all <- tibble()

for (m in metrics) {
  
  # Kruskal–Wallis
  kt <- kruskal.test(A[[m]] ~ A$Group)
  KW_all <- bind_rows(KW_all, tibble(
    Metric = m,
    KW_stat = as.numeric(kt$statistic),
    KW_p = kt$p.value
  ))
  
  # Dunn (always compute, like you wanted)
  dunn_res <- dunnTest(A[[m]] ~ A$Group, method="bh")$res
  dunn_res$Metric <- m
  DUNN_all <- bind_rows(DUNN_all, dunn_res)
}

write.csv(KW_all,  "Species_KruskalWallis_ALL.csv", row.names=FALSE)
write.csv(DUNN_all,"Species_DunnPostHoc_ALL.csv",  row.names=FALSE)

cat("Saved full Kruskal–Wallis and Dunn tables\n")

###############################################################################
### 7. MULTIVARIABLE ROBUST REGRESSION (HC3)
###############################################################################
A$Group <- factor(A$Group, levels=c("Healthy","Adenoma","Cancer"))
A$Sex   <- factor(A$Sex)

LM_out <- list()

for (m in metrics) {
  
  model <- lm(as.formula(paste0(m, " ~ Group + Age + BMI + Sex")), data=A)
  
  robust_se <- sqrt(diag(vcovHC(model, "HC3")))
  
  tbl <- tidy(model)
  tbl$robust_se <- robust_se
  tbl$robust_p  <- 2*(1 - pnorm(abs(tbl$estimate / robust_se)))
  tbl$Metric <- m
  
  LM_out[[m]] <- tbl
}

LM_df <- bind_rows(LM_out)
write.csv(LM_df, "Species_Alpha_LM_Robust.csv", row.names=FALSE)

cat("✔ Saved multivariable regression (HC3)\n")

###############################################################################
cat("\n Species Alpha Diversity analysis completed successfully.\n")
###############################################################################

