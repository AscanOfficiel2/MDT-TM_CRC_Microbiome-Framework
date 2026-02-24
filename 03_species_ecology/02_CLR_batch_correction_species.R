################################################################################
## Script Created BY Suleiman Aminu and AbdulAziz Ascandari
## FINAL PUBLICATION-READY VERSION
## Includes: metadata cleaning, TPM filtering, CLR, diagnostics,
## batch correction, and publication-quality PCA/UMAP/tSNE.
################################################################################

set.seed(43)

# --- Load libraries ---
library(limma)
library(compositions)
library(vegan)
library(ggplot2)
library(gridExtra)
library(umap)
library(Rtsne)
library(dplyr)
library(shadowtext)

################################################################################
# PUBLICATION-READY THEME + GROUP COLORS
################################################################################

pub_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.4),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

group_colors <- c(
  "Healthy" = "#ffcc00",
  "Adenoma" = "#984ea3",
  "Cancer"  = "#4daf4a"
)

################################################################################
# STEP 1 — LOAD TPM + METADATA AND CLEAN ALL NA
################################################################################

tpm  <- read.csv("TPM_matrix_arc.csv", row.names = 1, check.names = FALSE)
meta_raw <- read.csv("metadata_MAGs.csv", check.names = FALSE)

# Remove all NA-containing rows
meta_clean <- meta_raw[complete.cases(meta_raw), ]

# Align metadata to TPM columns
meta_clean <- meta_clean[match(colnames(tpm), meta_clean$Sample_ID), ]
meta_clean <- meta_clean[!is.na(meta_clean$Sample_ID), ]

write.csv(meta_clean, "metadata_final_matched.csv", row.names = FALSE)
cat("Saved metadata_final_matched.csv\n")

################################################################################
# STEP 2 — FILTER TPM TO MATCH CLEAN METADATA & SAVE RAW FILTERED TPM
################################################################################

tpm_filtered <- tpm[, meta_clean$Sample_ID, drop = FALSE]
tpm_filtered <- tpm_filtered[, match(meta_clean$Sample_ID, colnames(tpm_filtered))]

stopifnot(all(colnames(tpm_filtered) == meta_clean$Sample_ID))

tpm_filtered[is.na(tpm_filtered)] <- 0

write.csv(tpm_filtered, "TPM_raw_filtered.csv")
cat("Saved TPM_raw_filtered.csv\n")

################################################################################
# STEP 3 — CLR TRANSFORMATION
################################################################################

pseudocount <- 1e-6
tpm_clr <- t(apply(tpm_filtered + pseudocount, 2, clr))
#write.csv(tpm_clr, "TPM_clr.csv", quote = FALSE)

cat("Saved TPM_clr.csv\n")

################################################################################
# STEP 4 — ENSURE SAMPLES = ROWS
################################################################################

if (nrow(tpm_clr) == nrow(meta_clean)) {
  counts <- tpm_clr
} else {
  counts <- t(tpm_clr)
}

# Re-match metadata
meta_clean <- meta_clean[match(rownames(counts), meta_clean$Sample_ID), ]

################################################################################
# STEP 5 — REMOVE ZERO-VARIANCE TAXA
################################################################################

var_ok <- apply(counts, 2, function(x) var(x, na.rm = TRUE) > 1e-12)
counts_filt <- counts[, var_ok]

write.csv(counts_filt, "TPM_clr_filtered_zero_var.csv")
counts_t <- counts_filt

################################################################################
# --- PLOTTING FUNCTION FOR ALL POINT PLOTS ---
################################################################################

plot_points <- function(df, xvar, yvar, title) {
  ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(aes(fill = Group),
               color = "black", shape = 21,
               size = 3.5, stroke = 0.7, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    labs(title = title) +
    pub_theme
}

################################################################################
# STEP 6 — PCA BEFORE
################################################################################

pca <- prcomp(counts_t, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Group <- meta_clean$Group

p_before <- plot_points(pca_df, "PC1", "PC2", "PCA Before Batch Correction")

ggsave("PCA_Before_Correction.png", p_before,
       width = 5, height = 4, dpi = 600, bg = "white")

################################################################################
# STEP 7 — PERMANOVA BEFORE (EUCLIDEAN)
################################################################################

adonis_results <- list(
  Group     = adonis2(counts_t ~ Group,       data=meta_clean, method="euclidean"),
  Project   = adonis2(counts_t ~ Project,     data=meta_clean, method="euclidean"),
  Country   = adonis2(counts_t ~ Country,     data=meta_clean, method="euclidean"),
  Continent = adonis2(counts_t ~ Continent,   data=meta_clean, method="euclidean"),
  Center    = adonis2(counts_t ~ Center_Name, data=meta_clean, method="euclidean"),
  Instrument= adonis2(counts_t ~ Instrument,  data=meta_clean, method="euclidean")
)

adonis_df <- bind_rows(lapply(names(adonis_results), function(x){
  df <- as.data.frame(adonis_results[[x]])
  df$Factor <- x
  df
}))

permanova_clean <- adonis_df %>%
  filter(!is.na(F)) %>%         # keep only rows with valid F-statistic
  group_by(Factor) %>%
  slice(1) %>%                  # keep only the first row (the main test)
  ungroup()

# Save cleaned PERMANOVA table
write.csv(permanova_clean, "PERMANOVA_before_cleaned.csv", row.names = FALSE)
cat("Saved cleaned PERMANOVA table: PERMANOVA_before_cleaned.csv\n")


################################################################################
# STEP 8 — UMAP BEFORE
################################################################################

umap_res <- umap(counts_t)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df$Group <- meta_clean$Group

u_before <- plot_points(umap_df, "UMAP1", "UMAP2", "UMAP Before Batch Correction")
ggsave("UMAP_Before_Correction.png", u_before,
       width=5, height=5, dpi=600, bg="white")

################################################################################
# STEP 9 — t-SNE BEFORE
################################################################################

set.seed(42)
tsne_res <- Rtsne(counts_t, dims=2, perplexity=30)
tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1","tSNE2")
tsne_df$Group <- meta_clean$Group

ts_before <- plot_points(tsne_df, "tSNE1", "tSSE2", "t-SNE Before Batch Correction")

# Fix variable name issue
ts_before <- plot_points(tsne_df, "tSNE1", "tSNE2", "t-SNE Before Batch Correction")

ggsave("tSNE_Before.png", ts_before,
       width=5, height=5, dpi=600, bg="white")

################################################################################
# STEP 10 — BATCH CORRECTION USING LIMMA
################################################################################

design <- model.matrix(~ Group, data=meta_clean)
covars <- model.matrix(~ Instrument + Center_Name, data=meta_clean)[,-1]

# Version 1
counts_corrected_v1 <- removeBatchEffect(
  t(counts_t), covariates = covars, design = design)
counts_corrected_v1 <- t(counts_corrected_v1)
write.csv(counts_corrected_v1, "TPM_clr_batch_corrected_v1.csv")

# Version 2 (add Project)
counts_corrected_v2 <- removeBatchEffect(
  t(counts_t), batch = meta_clean$Project,
  covariates = covars, design = design)
counts_corrected_v2 <- t(counts_corrected_v2)
write.csv(counts_corrected_v2, "TPM_clr_batch_corrected_v2.csv")

################################################################################
# STEP 11 — DIAGNOSTICS AFTER CORRECTION
################################################################################

run_diagnostics <- function(corrected, meta, prefix) {
  
  corrected_t <- corrected
  
  # PCA
  pca <- prcomp(corrected_t, scale.=TRUE)
  pca_df <- as.data.frame(pca$x)
  pca_df$Group <- meta$Group
  
  p_after <- plot_points(pca_df, "PC1", "PC2",
                         paste("PCA After Correction", prefix))
  ggsave(paste0("PCA_after_", prefix, ".png"),
         p_after, width=6, height=5, dpi=600, bg="white")
  
  # PERMANOVA
  ad <- list(
    Group     = adonis2(corrected_t ~ Group,       data=meta, method="euclidean"),
    Instrument= adonis2(corrected_t ~ Instrument,  data=meta, method="euclidean"),
    Center    = adonis2(corrected_t ~ Center_Name, data=meta, method="euclidean"),
    Project   = adonis2(corrected_t ~ Project,     data=meta, method="euclidean")
  )
  
  ad_df <- bind_rows(lapply(names(ad), function(x){
    df <- as.data.frame(ad[[x]])
    df$Factor <- x
    df
  }))
  #write.csv(ad_df, paste0("PERMANOVA_after_", prefix, ".csv"), row.names=FALSE)
  
  
  # ---- CLEAN PERMANOVA TABLE: keep only main test rows ----
  permanova_clean <- ad_df %>%
    filter(!is.na(F)) %>%   # keep only rows with actual test statistics
    group_by(Factor) %>%
    slice(1) %>%            # keep FIRST row per factor (the main test)
    ungroup()
  
  # Save clean table
  write.csv(permanova_clean, paste0("PERMANOVA_after_cleaned_", prefix, ".csv"),
            row.names = FALSE)
  
  cat("Saved cleaned PERMANOVA table for", prefix, "\n")
  
  
  # UMAP
  umap_res <- umap(corrected_t)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1","UMAP2")
  umap_df$Group <- meta$Group
  
  u_after <- plot_points(umap_df, "UMAP1", "UMAP2",
                         paste("UMAP After Correction", prefix))
  ggsave(paste0("UMAP_after_", prefix, ".png"),
         u_after, width=6, height=5, dpi=600, bg="white")
  
  # t-SNE
  ts <- Rtsne(corrected_t, dims=2, perplexity=30)
  tsdf <- as.data.frame(ts$Y)
  colnames(tsdf) <- c("tSNE1","tSNE2")
  tsdf$Group <- meta$Group
  
  ts_after <- plot_points(tsdf, "tSNE1", "tSNE2",
                          paste("t-SNE After Correction", prefix))
  ggsave(paste0("tSNE_after_", prefix, ".png"),
         ts_after, width=6, height=5, dpi=600, bg="white")
}

# Run diagnostics
run_diagnostics(counts_corrected_v1, meta_clean, "v1")
run_diagnostics(counts_corrected_v2, meta_clean, "v2")

################################################################################
# END OF SCRIPT
################################################################################

