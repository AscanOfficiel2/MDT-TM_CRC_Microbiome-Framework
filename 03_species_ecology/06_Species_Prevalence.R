######################################################################
# MODULE 3 — PREVALENCE, CORE SPECIES, SHARED/UNIQUE SPECIES
######################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(UpSetR)
  library(pheatmap)
})

######################################################################
# 1. Prepare TPM matrix and metadata
######################################################################

# Reload aligned metadata (from Module 1)
metadata <- read.csv("CRC_metadata_aligned_bact.csv", check.names = FALSE)
rownames(metadata) <- metadata$Sample_ID

# TPM matrix
tpm <- read.csv("TPM_raw_filtered_bact.csv", row.names = 1, check.names = FALSE)

# Transpose if species are rows
if (nrow(tpm) > ncol(tpm)) tpm <- t(tpm)

tpm[is.na(tpm)] <- 0

# Ensure same sample order
tpm <- tpm[metadata$Sample_ID, ]

######################################################################
# 2. Presence–Absence Matrix
######################################################################

pa_matrix <- (tpm > 0) * 1   # convert TPM > 0 into presence/absence (1/0)

# Long format for prevalence analysis
pa_long <- pa_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(
    -Sample_ID,
    names_to = "Species",
    values_to = "Present"
  ) %>%
  left_join(metadata, by = "Sample_ID")

######################################################################
# 3. Prevalence by Group (Healthy, Adenoma, Cancer)
######################################################################

prevalence_env <- pa_long %>%
  group_by(Group, Species) %>%
  summarise(Prevalence = mean(Present) * 100, .groups = "drop")

write.csv(prevalence_env, "Species_Prevalence_ByGroup.csv", row.names = FALSE)

#######################################################################
# 4. Categorize species as Core / Secondary / Peripheral (Corrected)
######################################################################
prevalence_by_env <- prevalence_env %>%
  mutate(Category = case_when(
    Prevalence >= 80 ~ "Core",        # High prevalence
    Prevalence >= 50 ~ "Secondary",   # Medium prevalence
    TRUE ~ "Peripheral"               # Low prevalence
  ))

write.csv(
  prevalence_by_env,
  "Species_Prevalence_Group_Thresholds.csv",
  row.names = FALSE
)


######################################################################
# 5. Global Prevalence (Across ALL samples) — Corrected
######################################################################

global_prev <- colSums(pa_matrix) / nrow(pa_matrix) * 100

species_global <- data.frame(
  Species = names(global_prev),
  Prevalence = global_prev
) %>%
  mutate(Category = case_when(
    Prevalence >= 50 ~ "Core",
    Prevalence >= 30 ~ "Secondary",
    TRUE ~ "Peripheral"
  ))

write.csv(species_global, "Species_Global_Prevalence.csv", row.names = FALSE)

cat("Correct prevalence categories computed (Core ≥50, Secondary 20–50).\n")


write.csv(species_global, "Species_Global_Prevalence.csv", row.names = FALSE)

######################################################################
# 6. Global Prevalence Plot
######################################################################

summary_global <- species_global %>%
  group_by(Category) %>%
  summarise(n_species = n(), .groups = "drop") %>%
  mutate(Pct = round(100 * n_species / sum(n_species), 1))

tiff("Prevalence_Global.tiff", width = 4, height = 3.5, units = "in", res = 600)
ggplot(summary_global, aes(x = Category, y = n_species, fill = Category)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(Pct, "%")),
            vjust = -0.5, size = 2, fontface = "bold") +
  scale_fill_manual(values = c(
    "Core"       = "#1B9E77",   # Teal
    "Secondary"  = "#D95F02",   # Orange
    "Peripheral" = "#7570B3"    # Purple
  )) +

  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1)
  ) +
  labs(y = "Number of Species", x = "")
dev.off()

######################################################################
# 7. Prevalence per Group (Stacked Barplot)
######################################################################

summary_env <- prevalence_by_env %>%
  group_by(Group, Category) %>%
  summarise(n_species = n(), .groups = "drop")

tiff("Prevalence_ByGroup.tiff", width = 4.5, height = 3.5, units = "in", res = 600)
ggplot(summary_env, aes(x = Group, y = n_species, fill = Category)) +
  geom_col(color = "black", position = "stack") +
  scale_fill_manual(values = c(
    "Core"       = "#1B9E77",   # Teal
    "Secondary"  = "#D95F02",   # Orange
    "Peripheral" = "#7570B3"    # Purple
  ))  +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1)
  ) +
  labs(y = "Number of Species")
dev.off()

######################################################################
# 8. Heatmap — Top 50 Core Species per Group
######################################################################

######################################################################
# Heatmap — Top 50 Core Species per Group (fixed dplyr masking issue)
######################################################################

# 1. Select top core species
top_core <- prevalence_by_env %>%
  dplyr::filter(Category == "Core") %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(desc(Prevalence)) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::ungroup()

# 2. Build heatmap matrix safely
heatmap_matrix <- prevalence_by_env %>%
  dplyr::filter(Species %in% top_core$Species) %>%
  dplyr::select(Group, Species, Prevalence) %>%   # << EXPLICIT dplyr::select
  tidyr::pivot_wider(names_from = Group, values_from = Prevalence, values_fill = 0) %>%
  tibble::column_to_rownames("Species")

# 3. Plot
tiff("Heatmap_CoreSpecies.tiff", width = 4, height = 2.5, units = "in", res = 600)

pheatmap(
  as.matrix(heatmap_matrix),
  colorRampPalette(c("#FFFFD9","#41B6C4", "#081D58"))(100),
  #color = colorRampPalette(c("white","#6baed6","#08306b"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 5,
  fontsize_col = 7,
  border_color = NA
)

dev.off()

cat("Heatmap of top core species saved successfully.\n")


######################################################################
######################################################################
# STEP 7 — UpSet Plot (Real Species Names) with Correct Group Colors
######################################################################

library(UpSetR)

# Convert TRUE/FALSE to 1/0 for UpSetR
species_env_upset <- as.data.frame(species_env_matrix * 1)

# Your exact palette:
group_colors_upset <- c(
  "Healthy" = "#ffcc00",   # Yellow
  "Adenoma" = "#984ea3",   # Purple
  "Cancer"  = "#4daf4a"    # Green
)

# Output UpSet plot
tiff("UpSet_Species_ByGroup.tiff",
     width = 7, height = 4.5, units = "in", res = 600)

UpSetR::upset(
  species_env_upset,
  sets = colnames(species_env_upset),
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = group_colors_upset[colnames(species_env_upset)], # correct colors
  main.bar.color = "black",
  matrix.color = "black",
  point.size = 4,
  line.size = 1,
  text.scale = 1.5,
  number.angles = 0,
  mb.ratio = c(0.6, 0.4)
)

dev.off()

cat("\n UpSet plot saved as: UpSet_Species_ByGroup.tiff\n")

# 9. UpSet Plot — Shared Species Among Groups
######################################################################
######################################################################
# STEP 1 — Use REAL species names from TPM
######################################################################

# Ensure species names are correct
species_names <- colnames(tpm)

# Ensure metadata is aligned
stopifnot(all(rownames(tpm) == metadata$Sample_ID))

######################################################################
# STEP 2 — Build presence/absence matrix directly from TPM
######################################################################

pa_matrix <- (tpm > 0) * 1   # 1 = present, 0 = absent

######################################################################
# STEP 3 — Build a Species × Group presence matrix
######################################################################

species_env_matrix <- data.frame(
  Species = species_names,
  Healthy = colSums(pa_matrix[metadata$Group == "Healthy", ]) > 0,
  Adenoma = colSums(pa_matrix[metadata$Group == "Adenoma", ]) > 0,
  Cancer  = colSums(pa_matrix[metadata$Group == "Cancer", ]) > 0
)

rownames(species_env_matrix) <- species_env_matrix$Species
species_env_matrix$Species <- NULL   # now rownames = real species names

######################################################################
# STEP 4 — UNIQUE species per group
######################################################################

unique_species <- list()

for (grp in colnames(species_env_matrix)) {
  
  others <- setdiff(colnames(species_env_matrix), grp)
  
  sp <- rownames(species_env_matrix)[
    species_env_matrix[, grp] &
      rowSums(species_env_matrix[, others, drop = FALSE]) == 0
  ]
  
  unique_species[[grp]] <- sp
  
  cat("\n==============================\n")
  cat("UNIQUE SPECIES IN", grp, ":\n")
  print(sp)
  cat("==============================\n")
  
  write.csv(
    data.frame(Species = sp),
    paste0("Unique_Species_", grp, ".csv"),
    row.names = FALSE
  )
}

######################################################################
# STEP 5 — Shared by ALL groups
######################################################################

shared_all <- rownames(species_env_matrix)[
  rowSums(species_env_matrix) == ncol(species_env_matrix)
]

cat("\n==============================\n")
cat("SPECIES SHARED BY ALL GROUPS:\n")
print(shared_all)
cat("==============================\n")

write.csv(
  data.frame(Species = shared_all),
  "Shared_Species_AllGroups.csv",
  row.names = FALSE
)

######################################################################
# STEP 6 — Pairwise shared species
######################################################################

pairwise_shared <- list()
groups <- colnames(species_env_matrix)

for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    
    g1 <- groups[i]
    g2 <- groups[j]
    
    sp <- rownames(species_env_matrix)[
      species_env_matrix[, g1] &
        species_env_matrix[, g2] &
        rowSums(species_env_matrix[, setdiff(groups, c(g1, g2)), drop = FALSE]) == 0
    ]
    
    pairname <- paste0(g1, "_vs_", g2)
    pairwise_shared[[pairname]] <- sp
    
    cat("\n==============================\n")
    cat("SHARED SPECIES BETWEEN", g1, "AND", g2, ":\n")
    print(sp)
    cat("==============================\n")
    
    write.csv(
      data.frame(Species = sp),
      paste0("Shared_", pairname, ".csv"),
      row.names = FALSE
    )
  }
}

cat("\n DONE — Real Species Names Restored, Unique & Shared Species Correctly Computed.\n")

