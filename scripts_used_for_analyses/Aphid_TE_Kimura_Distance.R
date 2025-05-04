# Title: Kimura Distance Plot for Transposable Element Analysis
# Author: From Yann Bourgois, odified by Sergio Olvera
# Description: This script calculates the genome-normalized TE abundance across Kimura substitution levels 
#              and generates a stacked barplot showing TE composition over divergence levels.

# -----------------------------
# Load required libraries
# -----------------------------
library(reshape)      # for melting data
library(ggplot2)      # for plotting
library(viridis)      # for color scales
library(hrbrthemes)   # for modern ggplot2 themes
library(tidyverse)    # for data manipulation
library(gridExtra)
library(dplyr)
library(purrr)

# -----------------------------
# PARAMETERS
# -----------------------------
genome_size <- 317000000  # specify genome size in base pairs

# -----------------------------
# LOAD DATA
# -----------------------------
# Load first TE dataset
KimuraDistance <- read.table("file.txt", header = TRUE)

# -----------------------------
# PLOT TOTAL TE CONTENT
# -----------------------------
# Normalize by genome size and melt the dataframe
kd_melt <- melt(KimuraDistance, id = "Div")
kd_melt$norm <- kd_melt$value / genome_size * 100

# Define fill color (only one category assumed here)
colours <- c("aquamarine1")

# Plot bar chart of normalized TE abundance across Kimura distances
ggplot(kd_melt, aes(fill = variable, y = norm, x = Div)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = colours) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  labs(fill = "") 

# Save plot
ggsave("SMIS.png", width = 7.6, height = 6.6, dpi = 300)

# Total TE content
total_TE_content <- sum(KimuraDistance[, 2:ncol(KimuraDistance)]) / genome_size
print(total_TE_content)

# -----------------------------
# PROCESSING CLASSIFIED TEs
# -----------------------------
# Load categorized TE data
KimuraDistance <- read.table("SpeciesGenome.divsum3", header = TRUE)

# Collapse ambiguous TE types under "FalseAmbiguous"
col_true <- which(grepl("True", colnames(KimuraDistance)))
KimuraDistance$FalseAmbiguous <- rowSums(KimuraDistance[, col_true])

# Keep only "False" classifications
col_false <- which(grepl("False", colnames(KimuraDistance)))
KimuraDistance <- KimuraDistance[, c(1, col_false)]

# Standardize column names
colnames(KimuraDistance) <- gsub("TE_type\\.", "DNA/", colnames(KimuraDistance))
colnames(KimuraDistance) <- gsub("TE_Type", "DNA/", colnames(KimuraDistance))
colnames(KimuraDistance) <- gsub("TE_Type", "Gypsy", colnames(KimuraDistance))  # adjust as needed

# Optional: manually rename remaining ambiguous columns if needed
# colnames(KimuraDistance) <- c("Div", "DNA/TEType", "LTR/TEType") # example format

# Merge similar categories (e.g., unclassified subfamilies)
if ("Retrotransposon/" %in% colnames(KimuraDistance)) {
  KimuraDistance$`Retrotransposon/` <- rowSums(KimuraDistance[, grep("Retrotransposon", colnames(KimuraDistance))])
}
if ("DNA/" %in% colnames(KimuraDistance)) {
  KimuraDistance$`DNA/` <- rowSums(KimuraDistance[, grep("DNA_", colnames(KimuraDistance))])
}

# Keep only meaningful TE categories
final_cols <- grep("/", colnames(KimuraDistance), value = TRUE)
KimuraDistance <- KimuraDistance[, c("Div", final_cols)]

# -----------------------------
# PLOT CLASSIFIED TE COMPOSITION
# -----------------------------
kd_melt <- melt(KimuraDistance, id = "Div")
kd_melt$norm <- kd_melt$value / genome_size * 100

ggplot(kd_melt, aes(fill = variable, y = norm, x = Div)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  labs(fill = "TE Category")

# Total TE content again (optional)
total_TE_content_classified <- sum(KimuraDistance[, 2:ncol(KimuraDistance)]) / genome_size
print(total_TE_content_classified)