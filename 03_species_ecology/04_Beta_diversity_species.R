###############################################################################
# SPECIES-LEVEL BETA DIVERSITY ANALYSIS (PCoA + PERMANOVA + PERMDISP + ENVFIT + BIOENV)
###############################################################################

library(tidyverse)
library(vegan)
library(ggplot2)
library(plotly)
library(parallel)
library(tibble)

set.seed(42)

###############################################################################
# STEP 1 — LOAD CLR SPECIES MATRIX & METADATA
###############################################################################

cat("Loading CLR batch-corrected species matrix...\n")

clr_species <- read.csv("TPM_clr_batch_corrected_v2.csv",
                        row.names = 1, check.names = FALSE)

meta <- read.csv("CRC_metadata_aligned_bact.csv",
                 row.names = 1, check.names = FALSE)

# Align samples
common <- intersect(rownames(clr_species), rownames(meta))
clr_species <- clr_species[common, , drop=FALSE]
meta <- meta[common, , drop=FALSE]

cat("Samples aligned:", nrow(meta), "samples\n")

###############################################################################
# STEP 1B — REMOVE SAMPLES WITH MISSING METADATA
###############################################################################

cat("Filtering samples with missing (Group, Age, BMI, Sex)...\n")

meta <- meta %>% drop_na(Group, Age, BMI, Sex)
clr_species <- clr_species[rownames(meta), , drop=FALSE]

cat("Remaining samples:", nrow(meta), "\n\n")


###############################################################################
# STEP 3 — AITCHISON DISTANCE + PCoA (3 axes)
###############################################################################

dist_species <- dist(clr_species, method="euclidean")

cat("Running PCoA (k=3 axes)…\n")
PCOA <- cmdscale(dist_species, eig=TRUE, k=3)

PCOA_df <- as.data.frame(PCOA$points) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(meta %>% rownames_to_column("Sample_ID"), by="Sample_ID")

colnames(PCOA_df)[2:4] <- c("PCoA1", "PCoA2", "PCoA3")

# % variance explained
var_explained <- round(PCOA$eig[1:3] / sum(PCOA$eig) * 100, 2)
cat("Variance explained (%):\n")
print(var_explained)

###############################################################################
# STEP 4 — 2D PCoA PLOT (Aitchison)
###############################################################################

COLORS <- c("Healthy"="#ffcc00",
            "Adenoma"="#984ea3",
            "Cancer" ="#4daf4a")

p_pcoa2d <- ggplot(PCOA_df, aes(PCoA1, PCoA2, fill=Group)) +
  geom_point(shape=21, size=3, alpha=0.85, color="black") +
  scale_fill_manual(values=COLORS) +
  theme_classic(base_size=14) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1.1),
    legend.position="bottom"
  ) +
  labs(
    title="PCoA (Aitchison Distance) — Species",
    x=paste0("PCoA1 (", var_explained[1], "%)"),
    y=paste0("PCoA2 (", var_explained[2], "%)")
  )

ggsave("Species_PCoA2D_Aitchison.png",
       p_pcoa2d, width=4, height=5, dpi=600)

cat("Saved 2D PCoA plot\n\n")

#############################################################################
###############################################################################
# STEP 5 — STATIC 3D PCoA (NO plotly required)
###############################################################################
library(scatterplot3d)

png("Species_PCoA3D_Aitchison.png", width=2200, height=2000, res=300)

colors_vec <- c("Healthy"="#ffcc00",
                "Adenoma"="#984ea3",
                "Cancer" ="#4daf4a")

col_map <- colors_vec[PCOA_df$Group]

scatterplot3d(
  x = PCOA_df$PCoA1,
  y = PCOA_df$PCoA2,
  z = PCOA_df$PCoA3,
  color = col_map,
  pch = 19,
  xlab = paste0("PCoA1 (", var_explained[1], "%)"),
  ylab = paste0("PCoA2 (", var_explained[2], "%)"),
  zlab = paste0("PCoA3 (", var_explained[3], "%)"),
  main = "3D PCoA (Species | Aitchison Distance)"
)

legend("topright",
       legend = names(colors_vec),
       col = colors_vec,
       pch = 19)

dev.off()

cat("Static 3D PCoA saved: Species_PCoA3D_Aitchison.png\n")


###############################################################################
# STEP 6 — PERMANOVA (adonis2)
###############################################################################

meta_df <- meta %>% rownames_to_column("Sample_ID")

perm <- adonis2(dist_species ~ Group + Age + BMI + Sex,
                data = meta_df,
                permutations = 999,
                parallel = detectCores() - 1)

write.csv(as.data.frame(perm), "Species_PERMANOVA.csv")

cat("✔ PERMANOVA complete\n\n")
############################################################################
###############################################################################
# STEP 6B — PAIRWISE PERMANOVA (Group contrasts)
###############################################################################

pairwise_permanova <- function(distmat, groups, prefix) {
  
  groups <- as.factor(groups)
  lv <- levels(groups)
  combos <- combn(lv, 2, simplify = FALSE)
  out <- list()
  
  for (cmp in combos) {
    g1 <- cmp[1]
    g2 <- cmp[2]
    
    subset_idx <- groups %in% c(g1, g2)
    dist_sub <- as.dist(as.matrix(distmat)[subset_idx, subset_idx])
    meta_sub <- data.frame(Group = droplevels(groups[subset_idx]))
    
    pw <- adonis2(dist_sub ~ Group, data = meta_sub,
                  permutations = 999,
                  parallel = detectCores() - 1)
    
    out[[paste(g1, "_vs_", g2, sep="")]] <- tibble(
      Comparison = paste(g1, "vs", g2),
      R2  = pw$R2[1],
      F   = pw$F[1],
      p   = pw$`Pr(>F)`[1]
    )
  }
  
  result <- bind_rows(out)
  write.csv(result, paste0(prefix, "_Pairwise_PERMANOVA.csv"), row.names = FALSE)
  cat("Pairwise PERMANOVA saved as: ", paste0(prefix, "_Pairwise_PERMANOVA.csv"), "\n\n")
}

# Run for species
pairwise_permanova(dist_species, meta_df$Group, "Species")


###############################################################################
# STEP 7 — PERMDISP (Group homogeneity of dispersions)
###############################################################################

disp <- betadisper(dist_species, meta$Group)
write.csv(as.data.frame(anova(disp)), "Species_PERMDISP.csv")

cat("PERMDISP complete\n\n")
###############################################################################
###############################################################################
# PERMDISP Visualization — Distances to Group Centroid (Publication-Ready)
###############################################################################

disp <- betadisper(dist_species, meta$Group)
disp_df <- data.frame(
  Distance = disp$distances,
  Group = meta$Group
)

# Colors
COLORS <- c("Healthy"="#ffcc00",
            "Adenoma"="#984ea3",
            "Cancer" ="#4daf4a")

p_disp <- ggplot(disp_df, aes(Group, Distance, fill = Group)) +
  geom_jitter(width = 0.2, alpha = 0.25, size = 2, color="black") +
  geom_boxplot(width = 0.55, alpha = 0.85, color = "black", outlier.shape = NA) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "PERMDISP — Distance to Group Centroid",
    y = "Distance to Centroid",
    x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face="bold", size=16, hjust=0.5),
    axis.text.x = element_text(angle=45, hjust=1, face="bold"),
    axis.text.y = element_text(size=12),
    panel.grid = element_blank(),
    panel.border = element_rect(color="black", fill=NA, linewidth=1.2),
    legend.position = "none"
  )

ggsave("Species_PERMDISP_Boxplot.png",
       p_disp, width=5.2, height=5, dpi=600)

cat("PERMDISP visualization saved: Species_PERMDISP_Boxplot.png\n")
###########################################################################################################################
###############################################################################
# PERMDISP Visualization — Boxplot + Black Diamond + Dashed Line
###############################################################################

disp <- betadisper(dist_species, meta$Group)

disp_df <- data.frame(
  Distance = disp$distances,
  Group = factor(meta$Group, levels = c("Healthy","Adenoma","Cancer"))
)

COLORS <- c("Healthy"="#ffcc00",
            "Adenoma"="#984ea3",
            "Cancer" ="#4daf4a")

# Compute mean values per group
mean_df <- disp_df %>%
  group_by(Group) %>%
  summarise(mean_dist = mean(Distance))

p_disp_line <- ggplot(disp_df, aes(x = Group, y = Distance, fill = Group)) +
  
  geom_boxplot(
    width = 0.55,
    alpha = 0.85,
    color = "black",
    linewidth = 1.2,
    outlier.shape = NA
  ) +
  
  geom_line(
    data = mean_df,
    aes(x = Group, y = mean_dist, group = 1),
    color = "black",
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  
  # BLACK DIAMOND MEAN POINT
  geom_point(
    data = mean_df,
    aes(x = Group, y = mean_dist),
    shape = 23,       # diamond
    size = 3,
    fill = "white",   # ← FULL BLACK CENTER
    color = "black",  # border
    stroke = 1.0
  ) +
  
  scale_fill_manual(values = COLORS) +
  
  labs(
    title = "PERMDISP — Distance to Group Centroid",
    y = "Distance to Centroid",
    x = ""
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", size=18, hjust=0.5),
    axis.title.y = element_text(face="bold", size=16),
    axis.text.x = element_text(angle=45, hjust=1, face="bold", size=10, colour = "black"),
    axis.text.y = element_text(size=10, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color="black", fill=NA, linewidth=1.3),
    legend.position = "none"
  )

print(p_disp_line)

ggsave("Species_PERMDISP_Boxplot_Mean_BlackDiamond.png",
       p_disp_line, width = 4, height = 5, dpi = 600)

cat("PERMDISP boxplot with black diamond saved\n")


###############################################################################
# STEP 8 — ENVFIT (Age, BMI, Sex)
###############################################################################

env <- meta_df %>% select(Age, BMI, Sex)

cat("📊 Running ENVFIT…\n")

fit <- envfit(PCOA, env, permutations=999)

sink("Species_ENVFIT.txt"); print(fit); sink()

cat("ENVFIT saved\n\n")

###############################################################################
# STEP 9 — BIOENV (Best predictors of community structure)
###############################################################################

cat("Running BIOENV…\n")

env_num <- env %>%
  mutate(Sex_numeric = ifelse(Sex=="Male",1,0)) %>%
  select(Age, BMI, Sex_numeric)

bio <- bioenv(dist_species, env_num, method="spearman")

best_vars <- ifelse(is.null(bio$best),"None", paste(bio$best, collapse=" + "))
best_cor  <- ifelse(is.null(bio$bestcor), NA, bio$bestcor)

bio_df <- tibble(
  Best_Model = best_vars,
  Correlation = best_cor
)

write.csv(bio_df, "Species_BIOENV.csv", row.names=FALSE)

cat("✔ BIOENV complete\n\n")

###############################################################################
cat("\n Species Beta Diversity Pipeline Completed Successfully!\n")
###############################################################################
