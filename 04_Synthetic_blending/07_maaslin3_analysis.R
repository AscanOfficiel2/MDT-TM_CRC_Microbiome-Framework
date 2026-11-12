library(maaslin3)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ------------------------------------------------------------
# FILE PATHS
# ------------------------------------------------------------
species_tpm_file <- "TPM_raw_filtered_bact.csv"
metadata_file    <- "CRC_metadata_with_MDT_state.csv"
output_dir       <- "Maaslin3_MDT_state"

dir.create(output_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# LOAD SPECIES TPM MATRIX
# (rows = species, cols = samples)
# ------------------------------------------------------------
species_tpm_raw <- read_csv(
  species_tpm_file,
  show_col_types = FALSE
) %>% as.data.frame()

rownames(species_tpm_raw) <- species_tpm_raw[[1]]
species_tpm_raw <- species_tpm_raw[, -1]

# ------------------------------------------------------------
# TRANSPOSE → samples as rows
# ------------------------------------------------------------
species_tpm <- as.data.frame(t(species_tpm_raw))
species_tpm$Sample_ID <- rownames(species_tpm)
rownames(species_tpm) <- species_tpm$Sample_ID
species_tpm$Sample_ID <- NULL

cat("✔ Species TPM matrix shape:", dim(species_tpm), "\n")

# ------------------------------------------------------------
# LOAD METADATA
# ------------------------------------------------------------
meta <- read_csv(metadata_file, show_col_types = FALSE) %>%
  as.data.frame()

rownames(meta) <- meta$Sample_ID

# ------------------------------------------------------------
# RESTRICT TO ADENOMA + CANCER
# ------------------------------------------------------------
meta <- meta %>%
  filter(Group %in% c("Adenoma", "Cancer"))

species_tpm <- species_tpm[rownames(meta), , drop = FALSE]

cat("Samples used:", nrow(meta), "\n")

# ------------------------------------------------------------
# FACTOR SETUP
# ------------------------------------------------------------
meta$MDT_state <- factor(meta$MDT_state,
                         levels = c("pre_MDT", "post_MDT"))

# ------------------------------------------------------------
# RUN MAASLIN3
# ------------------------------------------------------------
fit <- maaslin3(
  input_data     = species_tpm,
  input_metadata = meta,
  output         = output_dir,
  
  formula        = "~ MDT_state",
  
  normalization  = "TSS",
  transform      = "LOG",
  correction     = "BH",
  
  standardize    = FALSE,
  augment        = TRUE,
  
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = FALSE,
  
  max_significance = 0.1,
  cores = 1
)

cat("MaAsLin3 MDT_state analysis completed successfully\n")
################################################################################
# ------------------------------------------------------------
# Prepare volcano data (MATCH MaAsLin summary logic)
# ------------------------------------------------------------
maaslin <- read_tsv(
  "Maaslin3_MDT_state/all_results.tsv",
  show_col_types = FALSE
)

# Keep abundance model only
maaslin <- maaslin %>%
  filter(model == "abundance") %>%
  mutate(
    neg_log10_q = -log10(qval_joint),
    significant = qval_joint <= 0.05
  )

# ------------------------------------------------------------
# Save ALL joint-significant species (abundance model)
# ------------------------------------------------------------
sig_all <- maaslin %>%
  filter(significant) %>%
  arrange(desc(neg_log10_q))

write_csv(
  sig_all,
  "Maaslin3_MDT_joint_significant_abundance_species.csv"
)

cat("✔ CSV saved: Maaslin3_MDT_joint_significant_abundance_species.csv\n")

# ------------------------------------------------------------
# Select top 4 from each side for labeling
# ------------------------------------------------------------
sig_df <- maaslin %>%
  filter(significant)

top_pos <- sig_df %>%
  filter(coef > 0) %>%
  slice_max(
    order_by = abs(coef) * neg_log10_q,
    n = 4,
    with_ties = FALSE
  )

top_neg <- sig_df %>%
  filter(coef < 0) %>%
  slice_max(
    order_by = abs(coef) * neg_log10_q,
    n = 4,
    with_ties = FALSE
  )

label_df <- bind_rows(top_pos, top_neg)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(maaslin,
            aes(x = coef, y = neg_log10_q)) +
  
  geom_point(aes(color = significant),
             size = 3,
             alpha = 0.6) +
  
  scale_color_manual(values = c("grey70", "#2C7BB6")) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_text_repel(
    data = label_df,
    aes(label = feature),
    size = 2.3,          # slightly smaller font
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  
  labs(
    title = "Differentially Abundant Species — MDT boundary",
    subtitle = "Adenoma → Cancer transition (joint FDR)",
    x = "Effect Size (log fold-change)",
    y = expression(-log[10]("Joint FDR q-value"))
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    legend.position = "none",
    
    # Remove default half frame
    axis.line = element_blank(),
    
    # Full thick border
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.8
    ),
    
    # Thicker ticks
    axis.ticks = element_line(linewidth = 1.2),
    
    # Slightly stronger axis text
    axis.text = element_text(size = 11),
    
    # Bold axis titles
    axis.title = element_text(face = "bold")
  )

ggsave("Maaslin3_MDT_Volcano_Annotated_joint.png",
       p,
       width = 5,
       height = 5,
       dpi = 600)

print(p)
#########################################################################
###Healthy to Adenoma
