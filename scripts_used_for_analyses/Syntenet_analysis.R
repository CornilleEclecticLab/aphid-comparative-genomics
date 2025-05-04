---
  title: "Inference and analysis of synteny networks for aphid genomes"
author: "Fabricio Almeida-Silva modified by Aurelie Mesnil"
---
  
  # Load required packages
  library(here)
library(syntenet)
library(tidyverse)
library(readxl)
library(patchwork)

# Load GR and OR gene lists from Excel file
gr_path <- here("data", "Gene_list.xlsx")
gr_list <- excel_sheets(gr_path) %>% set_names() %>% map(read_excel, path = gr_path)

or_path <- here("data", "Gene_list.xlsx")
or_list <- excel_sheets(or_path) %>% set_names() %>% map(read_excel, path = or_path)

# Load GFF annotations and protein FASTA sequences
annot_files <- dir(here("data", "annotation"), full.names = TRUE)
annot <- lapply(annot_files, function(x) {
  gr <- rtracklayer::import(x, feature.type = "gene")
  gr[, c("source", "type", "ID", "Name")]
}) %>% setNames(gsub("_.*", "", basename(annot_files)))
annot <- GenomicRanges::GRangesList(annot)

seq <- fasta2AAStringSetlist(here("data", "sequences"))
names(seq) <- gsub("_.*", "", names(seq))

# Process input for syntenet
pdata <- process_input(seq, annot)
pdata$seq <- lapply(pdata$seq, function(x) x[width(x) > 15])
save(pdata, compress = "xz", file = here("products", "result_files", "pdata.rda"))

# DIAMOND similarity search and synteny inference
load(here("products", "result_files", "pdata.rda"))
diamond <- run_diamond(seq = pdata$seq)
aphid_syntenet <- infer_syntenet(blast_list = diamond, annotation = pdata$annotation, outdir = file.path(tempdir(), "aphid_syn_output"))
save(aphid_syntenet, compress = "xz", file = here("products", "result_files", "aphid_syntenet.rda"))

# Cluster synteny network
clusters <- cluster_network(aphid_syntenet)
save(clusters, compress = "xz", file = here("products", "result_files", "clusters.rda"))

# Generate phylogenomic profiles
profiles <- phylogenomic_profile(clusters)

# Load species metadata and prepare order
tree_order <- c("MCE", "MPE", "DPL", "DNO", "API", "SMI", "PNI", "AGL", "AGO", "RMA", "RPA", "ELA", "DVI")
species_metadata <- read_excel(here("data", "species_metadata.xlsx")) %>%
  janitor::clean_names() %>%
  filter(!is.na(acronyms)) %>%
  mutate(acronyms = str_sub(acronyms, 1, 3)) %>%
  select(species, acronyms, origin, life_cycle, reproduction, host_specialization) %>%
  inner_join(data.frame(acronyms = tree_order), by = "acronyms")

species_order <- setNames(species_metadata$acronyms, species_metadata$species)

# Identify syntenic OR and GR genes
gr_genes <- unique(Reduce(rbind, gr_list)$`Gene ID`)
or_genes <- unique(Reduce(rbind, or_list)$`Gene ID`)

gr_clusters <- clusters %>% mutate(Gene = str_replace_all(Gene, "[a-zA-Z]{3,5}_", "")) %>% filter(Gene %in% gr_genes)
or_clusters <- clusters %>% mutate(Gene = str_replace_all(Gene, "[a-zA-Z]{3,5}_", "")) %>% filter(Gene %in% or_genes)

gr_or_clusters <- inner_join(rename(gr_clusters, Gene_GR = Gene), rename(or_clusters, Gene_OR = Gene))

chemosensory_syntenic <- data.frame(
  GR_syntenic = nrow(gr_clusters),
  GR_syntenic_percentage = nrow(gr_clusters) / length(gr_genes),
  OR_syntenic = nrow(or_clusters),
  OR_syntenic_percentage = nrow(or_clusters) / length(or_genes)
)

# Plot phylogenomic profiles
profile_gr <- profiles[unique(gr_clusters$Cluster), ]
profile_or <- profiles[unique(or_clusters$Cluster), ]

p_prof_gr <- plot_profiles(profile_gr, species_annotation = species_metadata %>% select(acronyms, specialization = host_specialization), cluster_species = species_order, show_colnames = TRUE, main = "Phylogenomic profiles of GR-encoding genes")

p_prof_or <- plot_profiles(profile_or, species_annotation = species_metadata %>% select(acronyms, specialization = host_specialization), cluster_species = species_order, show_colnames = TRUE, main = "Phylogenomic profiles of OR-encoding genes")

# Optional: plot by origin
p_prof_gr_origin <- plot_profiles(profile_gr, species_annotation = species_metadata %>% select(acronyms, origin), cluster_species = species_order, show_colnames = TRUE, main = "Phylogenomic profiles of GR-encoding genes", border_color = "black")

# Synteny network visualizations
genes_df <- unique(c(aphid_syntenet$Anchor1, aphid_syntenet$Anchor2)) %>% as.data.frame() %>% rename(gene = 1)
color_by <- inner_join(gene_df <- genes_df %>% mutate(acronyms = str_replace_all(gene, "_.*", "")), species_metadata)

c_gr <- filter(clusters, Cluster %in% unique(gr_clusters$Cluster))
gr_netplot <- plot_network(aphid_syntenet, clusters = c_gr, cluster_id = unique(c_gr$Cluster), color_by = color_by[, c("gene", "species")])

c_or <- filter(clusters, Cluster %in% unique(or_clusters$Cluster))
or_netplot <- plot_network(aphid_syntenet, clusters = c_or, cluster_id = unique(c_or$Cluster), color_by = color_by[, c("gene", "species")])

# Combine GR and OR network plots
combined_netplot <- wrap_plots(
  gr_netplot + ggtitle("Syntenic relationships among GR-encoding genes"),
  or_netplot + ggtitle("Syntenic relationships among OR-encoding genes")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_netplot
