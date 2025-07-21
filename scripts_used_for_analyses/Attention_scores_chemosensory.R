---
  title: "Attention scores analysis for chemosensory genes"
author: "Aurelie Mesnil"
---

  # Load required packages

library(rhdf5)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(here)
library(tidyverse)

rm(list = ls())


# ---- Extract attention scores from h5 files ----

# Path of the HDF5 file
root_dir <- "results/attentions_chemosensory_gene" 
h5_files <- list.files(path = root_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)

list.files(path = root_dir)

# Extraction function 
extract_attention_from_file <- function(h5_file) {
  gene_names <- h5ls(h5_file)$name
  file_id <- basename(h5_file)  # Ou extraire un nom plus pertinent si besoin
  
  attention_list <- lapply(gene_names, function(gene) {
    scores <- h5read(h5_file, gene)
    data.frame(
      Gene = gene,
      Position = seq_along(scores),
      Attention = scores,
      File = file_id
    )
  })
  
  do.call(rbind, attention_list)
}

all_attention_data <- lapply(h5_files, extract_attention_from_file)
attention_df <- bind_rows(all_attention_data) %>% 
  mutate(Gene_family = sub("_.*", "", File)) %>% 
  mutate(Chemotype = ifelse(grepl("OR", Gene_family), "OR", "GR")) 

unique(attention_df$Gene_family) %>% length

attention_df %>% ggplot(aes(x = Attention))+
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~ Chemotype)+
  geom_vline(xintercept = 0.5, colour = "red")+
  theme_minimal()



# ---- Add selection signal information on each site ---

list_files = list.files(path = "results", pattern = "Real_positions_")
Positives_sites = do.call(rbind, lapply(seq_along(list_files), function(i){
  gene_family = gsub("Real_positions_", "", list_files[i])
  gene_family = gsub(".csv", "", gene_family)
  dt = fread(here("results",list_files[i]), header = T) %>%
    pivot_longer(-1, names_to = "alignment_pos") %>% 
    mutate(gene_family = gene_family) %>% 
    dplyr::rename(Gene_id = 1)
  
})) %>% mutate(positive_keys = paste(value, Gene_id, sep = "_"))

attention_df = attention_df %>% mutate(key = paste(Position, Gene, sep = "_"))
attention_df = attention_df %>% 
  mutate(type_site = ifelse(key %in% Positives_sites$positive_keys, "Positive selection", "Neutral"))

count_gene_by_fam = attention_df %>% group_by(Gene_family, Gene) %>% count() 

GR06 = attention_df %>% filter(Gene_family == "GR06")
  
# ---- Add association between gene and TE information ---

ChemTE_association = fread("results/Enriched_te.csv") %>% distinct(Gene_id) %>% pull(Gene_id)

attention_df = attention_df %>% mutate(TE_association = ifelse(Gene %in% ChemTE_association, "TE associated", 'Not TE associated'))


# ---- Count by group ----

attention_df %>%
  group_by(Chemotype, TE_association, type_site) %>%
  count() %>% print()


# ---- Preparation of comparisons and statistical tests ----
attention_df$type_site <- factor(attention_df$type_site, levels = c("Neutral", "Positive selection"))

# ---- Make plots ----
plot_data <- attention_df %>%
  group_by(Chemotype, TE_association, type_site) %>%
  mutate(N = n()) %>%
  ungroup()

plot_data = filter(plot_data, Gene_family %in% c("GR03", "GR04", "GR06", "OR20"))

p0 = ggplot(plot_data, aes(x = type_site, y = Attention))+
  geom_boxplot(outlier.shape = NA, fill = "lightgrey", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  stat_summary(fun = median, geom = "point", shape = 21, size = 3, fill = "red") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3,
                     comparisons = list(c("Neutral", "Positive selection")),
                     label.y = 1.1) +
  stat_summary(fun.data = function(x) {
    data.frame(y = 1.05, label = paste("n =", length(x)))
  }, geom = "text", size = 3) +
  # facet_grid(~ Chemotype)+
  facet_grid(~ Gene_family)+
  labs(x = "", y = "Attention score") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p0


p1 = ggplot(plot_data, aes(x = type_site, y = Attention)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgrey", color = "black") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 3, fill = "red") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3,
                     comparisons = list(c("Neutral", "Positive selection")),
                     label.y = 1.1) +
  stat_summary(fun.data = function(x) {
    data.frame(y = 1.05, label = paste("n =", length(x)))
  }, geom = "text", size = 3) +
  facet_grid( ~ Chemotype) +
  labs(x = "", y = "Attention score") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

p1



p2 = attention_df %>%
  ggplot(aes(x = TE_association, y = Attention)) +
  geom_boxplot(fill = "lightgrey", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "red") +
  
  stat_compare_means(
    method = "wilcox.test", label = "p.format", size = 3,
    comparisons = list(c("TE associated", 'Not TE associated')),
    label.y = 1.1, )+
  
  stat_summary(fun.data = function(x) {
    data.frame(y = 1.05, label = paste("n =", length(x)))
  }, geom = "text", size = 3) +
  facet_wrap(~ Chemotype, nrow = 1) +
  labs(
    x = "",
    y = "Attention score") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

p2

p3 = attention_df %>%
  filter(type_site == "Positive selection") %>% 
  ggplot(aes(x = TE_association, y = Attention)) +
  geom_boxplot(fill = "lightgrey", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "red") +
  
  stat_compare_means(
    method = "wilcox.test", label = "p.format", size = 3,
    comparisons = list(c("TE associated", 'Not TE associated')),
    label.y = 1.1, )+
  
  stat_summary(fun.data = function(x) {
    data.frame(y = 1.05, label = paste("n =", length(x)))
  }, geom = "text", size = 3) +
  facet_wrap(~ Chemotype, nrow = 1) +
  labs(
    x = "",
    y = "Attention score") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

p3


# ---- Fisher test to test enrichment of positive selection sites in TE associated genes ----

attention_df %>%
  group_by(Chemotype) %>%
  group_map(~ fisher.test(table(.x$TE_association, .x$type_site)))



data = attention_df %>%
  filter(Chemotype == "OR")

ct = table(data$TE_association, data$type_site)
# Valeurs attendues manuellement
expected <- outer(rowSums(ct), colSums(ct)) / sum(ct)

expected

FT = fisher.test(ct)
