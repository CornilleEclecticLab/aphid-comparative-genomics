---
  title: "Study association between enriched TE and chemosensory genes"
author: "Adaptation of TEgrip https://github.com/marieBvr/TEs_genes_relationship_pipeline"
---

# Load required libraries
library(tidyverse)
library(data.table)
library(seqinr)
library(reshape2)
library(plotly)
library(ggh4x)
library(RColorBrewer)

# Define working path and input data
path <- "/home/aurelie/Dropbox/Post-Doc_ECLECTIC/8_Aphids_comparative_genomic/tegrip_aphids/"

# List of TEs enriched according to LOLA
list_TE_lola <- c(
  "AGLY_TEdenovoGr-B-G3101-Map3", "DPLA_TEdenovoGr-B-G149-Map20_reversed", "DPLA_TEdenovoGr-B-G161-Map14",
  "DPLA_TEdenovoGr-B-G4075-Map8_reversed", "DPLA_TEdenovoGr-B-G4329-Map4", "DPLA_TEdenovoGr-B-G77-Map20",
  "DVIT_TEdenovoGr-B-G6314-Map6", "DVIT_TEdenovoGr-B-G3298-Map9_reversed", "DVIT_TEdenovoGr-B-G4429-Map3",
  "DVIT_TEdenovoGr-B-G678-Map3", "DVIT_TEdenovoGr-B-G7829-Map4", "DVIT_TEdenovoGr-B-G9099-Map3",
  "MCER_TEdenovoGr-B-G1335-Map3", "RPAD_TEdenovoGr-B-G1099-Map3", "RPAD_TEdenovoGr-B-G1773-Map4"
)

# List of species used in the study
sp_list <- c("AGLY", "AGOS", "DNOX", "DPLA", "DVIT", "ELAN", "MCER", "MPER", "PNIG", "RMAI", "RPAD", "SMIS")

# Define a custom color palette
colors_list <- brewer.pal(n = 6, name = "Spectral")

# TEgrip-like function: for each TE-gene pair, calculate spatial relationship
calcul_distance <- function(TE, Gene_family) {
  res <- data.frame(
    Chr = character(), TE_id = character(), TE_strand = character(), TE_start = numeric(), TE_end = numeric(),
    Gene_id = character(), Gene_strand = character(), Gene_start = numeric(), Gene_end = numeric(),
    distance1 = numeric(), distance2 = numeric(), distance3 = numeric(), distance4 = numeric(),
    relationship = character(), Gene_fam = character(), stringsAsFactors = FALSE
  )
  
  species <- gsub("_TEdenovoGr.*", "", TE)
  
  # Load TE coordinates
  TE_input <- fread(paste0(path, "data/", species, "_TE_annotation.tsv")) %>%
    select(seqid, start, end, strand, class, TEname) %>%
    filter(TEname == TE)
  
  # Load gene coordinates
  gene_input <- fread(paste0(path, "data/", species, "_", Gene_family, "gene_annotation.tsv")) %>%
    select(chromosome, start, end, strand, ID)
  
  liste_genes <- unique(gene_input$ID)
  
  # Iterate through each gene to find nearby TEs and determine spatial relationship
  for (gene in liste_genes) {
    gene_info <- filter(gene_input, ID == gene)
    chromosome <- gene_info$chromosome
    start <- gene_info$start
    end <- gene_info$end
    strand <- gene_info$strand
    
    TE_in_same_chr <- filter(TE_input, seqid == chromosome)
    
    if (nrow(TE_in_same_chr) > 0) {
      tmp <- data.frame(
        Chr = chromosome,
        TE_id = TE,
        TE_strand = TE_in_same_chr$strand,
        TE_start = TE_in_same_chr$start,
        TE_end = TE_in_same_chr$end,
        Gene_start = start,
        Gene_end = end,
        Gene_strand = strand,
        distance1 = TE_in_same_chr$start - end,
        distance2 = start - TE_in_same_chr$end,
        distance3 = end - TE_in_same_chr$end,
        distance4 = start - TE_in_same_chr$start,
        Gene_id = gene,
        relationship = "",
        Gene_fam = Gene_family
      )
      
      # Classify spatial relationship between TE and gene
      tmp <- tmp %>%
        rowwise() %>%
        mutate(relationship = case_when(
          distance1 < 0 & distance2 < 0 & distance3 >= 0 & distance4 <= 0 ~ "subset",
          distance1 < 0 & distance2 < 0 & distance3 <= 0 & distance4 >= 0 ~ "superset",
          distance1 < 0 & distance2 < 0 & distance3 < 0 & distance4 < 0 ~ "downstream_overlap",
          distance1 < 0 & distance2 < 0 & distance3 > 0 & distance4 > 0 ~ "upstream_overlap",
          distance1 > 0 & distance2 < 0 & distance3 < 0 & distance4 < 0 ~ "downstream",
          distance1 < 0 & distance2 > 0 & distance3 > 0 & distance4 > 0 ~ "upstream",
          TRUE ~ ""
        )) %>%
        ungroup() %>%
        mutate(relationship = ifelse(
          (relationship %in% c("upstream", "downstream")) &
            pmin(abs(distance1), abs(distance2), abs(distance3), abs(distance4)) > 10000,
          "", relationship
        )) %>%
        filter(relationship != "")
      
      if (nrow(tmp) > 0) {
        res <- bind_rows(res, tmp)
      }
    }
  }
  return(res)
}

# Process all TE and gene families
results_list <- list()

for (TE in list_TE_lola) {
  res_OR <- calcul_distance(TE, "OR") %>%
    filter(relationship != "") %>%
    select(-starts_with("distance"))
  
  results_list[[length(results_list) + 1]] <- res_OR
  
  res_GR <- calcul_distance(TE, "GR") %>%
    filter(relationship != "") %>%
    select(-starts_with("distance"))
  
  results_list[[length(results_list) + 1]] <- res_GR
}

# Combine all results
dt_all <- rbindlist(results_list)
write.table(dt_all, file = paste0(path, "results/Enriched_te.csv"), row.names = FALSE, sep = "\t")

# Count how many genes are associated to each TE
Count_genes <- dt_all %>%
  filter(TE_id %in% list_TE_lola) %>%
  mutate(species = str_remove_all(TE_id, "_TEdenovoGr.*")) %>%
  group_by(TE_id, Gene_fam, relationship, species) %>%
  count()

write.table(Count_genes, file = paste0(path, "results/Count_genes.csv"), row.names = FALSE, sep = "\t")

# Count how many TEs are associated to each gene
Count_TE <- dt_all %>%
  filter(TE_id %in% list_TE_lola) %>%
  mutate(species = str_remove_all(TE_id, "_TEdenovoGr.*")) %>%
  group_by(Gene_id, Gene_fam, relationship, species) %>%
  count()

write.table(Count_TE, file = paste0(path, "results/Count_TE.csv"), row.names = FALSE, sep = "\t")

# Create figure showing gene counts per TE and relationship category
fig3_A <- Count_genes %>%
  ggplot(aes(x = n, y = TE_id, fill = relationship)) +
  geom_bar(stat = "identity") +
  facet_grid(species ~ Gene_fam, scales = "free_y") +
  force_panelsizes(rows = c(1, 5, 6, 1, 2)) +
  scale_fill_manual(values = setNames(colors_list, 
    c("downstream", "downstream_overlap", "subset", "superset", "upstream", "upstream_overlap"))) +
  scale_y_discrete(labels = c(
    "AGLY_TEdenovoGr-B-G3101-Map3" = "AGLY-G3101-Map3 | TIR hAT",
    "DPLA_TEdenovoGr-B-G149-Map20_reversed" = "DPLA-G149-Map20 | TIR Tc1-Mariner",
    "DPLA_TEdenovoGr-B-G161-Map14" = "DPLA-G161-Map14 | class II",
    "DPLA_TEdenovoGr-B-G4075-Map8_reversed" = "DPLA-G4075-Map8 | SINE",
    "DPLA_TEdenovoGr-B-G4329-Map4" = "DPLA-G4329-Map4 | TIR hAT",
    "DPLA_TEdenovoGr-B-G77-Map20" = "DPLA-G77-Map20 | class II",
    "DVIT_TEdenovoGr-B-G6314-Map6" = "DVIT-G6314-Map6 | class II",
    "DVIT_TEdenovoGr-B-G3298-Map9_reversed" = "DVIT-G3298-Map9 | TIR hAT",
    "DVIT_TEdenovoGr-B-G4429-Map3" = "DVIT-G4429-Map3 | TIR hAT",
    "DVIT_TEdenovoGr-B-G678-Map3" = "DVIT-G678-Map3 | TIR Mutator",
    "DVIT_TEdenovoGr-B-G7829-Map4" = "DVIT-G7829-Map4 | TIR hAT",
    "DVIT_TEdenovoGr-B-G9099-Map3" = "DVIT-G9099-Map3 | TIR hAT",
    "MCER_TEdenovoGr-B-G1335-Map3" = "MCER-G1335-Map3 | SINE",
    "RPAD_TEdenovoGr-B-G1099-Map3" = "RPAD-G1099-Map3 | TIR Tc1-Mariner",
    "RPAD_TEdenovoGr-B-G1773-Map4" = "RPAD-G1773-Map4 | unclassified"
  )) +
  labs(x = "Gene count", y = "Enriched TEs") +
  theme_linedraw() +
  theme(strip.background = element_rect(fill = "grey80", colour = "grey30"),
        strip.text = element_text(color = "black"))

ggsave(fig3_A, filename = "results/fig3_A.png", width = 15, height = 30, units = "cm")

# Compute how many OR/GR genes are enriched in TEs
TEenriched_genes <- unique(dt_all$Gene_id)

Chemgenes_annotation <- do.call(rbind, lapply(sp_list, function(sp) {
  ORgenes <- fread(paste0(path, "data/", sp, "_ORgene_annotation.tsv")) %>%
    mutate(species = sp, gene_class = "OR")
  GRgenes <- fread(paste0(path, "data/", sp, "_GRgene_annotation.tsv")) %>%
    mutate(species = sp, gene_class = "GR")
  bind_rows(ORgenes, GRgenes)
})) %>%
  mutate(TE_enriched_gene = ifelse(ID %in% TEenriched_genes, "yes", "no"))

# Summarize TE association by species and gene class
Chemgenes_annotation %>%
  group_by(species, TE_enriched_gene, gene_class) %>%
  tally() %>%
  print(n = 32)