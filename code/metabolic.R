library(tidyverse)
library(ggh4x)

# ----------------------------
# File paths
# ----------------------------
annotation_fp <- "./data/metabolic_g/METABOLIC_result.xlsx/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv"
color_fp <- "/projects/p32449/color_palettes/color_dictionary_draft07.csv"
classification_fp <- "./data/gtdbtk_all/gtdbtk.bac120.summary.tsv"
sample_key_fp <- "./data/sharepoint_data/submission_key.csv"
ecofold_fp <- "./data/EcoFoldDB/annotated"
ecofold_colors_fp <- "/projects/p32449/color_palettes/EcoFoldDB_subcat_colors_df.rds"

# ----------------------------
# Read data
# ----------------------------
annotation_raw <- read_tsv(annotation_fp)
color_raw <- read_csv(color_fp)
classification_raw <- read_tsv(classification_fp)
sample_key <- read_csv(sample_key_fp)
ef_color_raw <- readRDS(ecofold_colors_fp) # data frame: Category, Subcat_collapsed, color, order
eco_files <- eco_files <- list.files(
  ecofold_fp,
  pattern = "_annotations\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# ----------------------------
# Read and merge EcoFoldDB outputs
# ----------------------------
eco_df <- eco_files %>%
  map_df(~ {
    fname <- basename(.x)
    sample <- sub("\\..*", "", fname)
    read_tsv(.x) %>% mutate(sample = sample, filename = fname)
  }) %>%
  left_join(sample_key, by = "sample") %>%
  mutate(subcat = case_when(
    str_detect(`Sub-category`, "Polyphenol") ~ "Polyphenol cycling",
    TRUE ~ `Sub-category`
  )) %>%
  left_join(classification_raw, by = c("sample" = "user_genome")) %>%
  separate(
    classification, 
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";", fill = "right", remove = FALSE
  ) %>%
  mutate(across(domain:species, ~ sub(".*__", "", .))) %>%
  left_join(ef_color_raw, by = c("subcat" = "Subcategory", "Category" = "Category"))%>%
  mutate(
    subcat = factor(subcat, levels = ef_color_raw$Subcategory[order(ef_color_raw$order)]),
    Category = factor(Category, levels = unique(ef_color_raw$Category[order(ef_color_raw$order)]))
  )

# ----------------------------
# Tidy METABOLIC annotation data
# ----------------------------
annotation_tidy <- annotation_raw %>%
  mutate(across(-c(Category:`Hmm detecting threshold`), as.character)) %>%
  pivot_longer(
    cols = -c(Category:`Hmm detecting threshold`),
    names_to = "bin_measure",
    values_to = "value"
  ) %>%
  separate(bin_measure, into = c("sample_bin", "measurement"), sep = " ", extra = "merge") %>%
  pivot_wider(names_from = "measurement", values_from = "value")

# ----------------------------
# Add taxonomy to METABOLIC data
# ----------------------------
annotation_classified <- annotation_tidy %>%
  mutate(sample = str_remove_all(sample_bin, "_scaffolds")) %>%
  left_join(classification_raw, by = c("sample" = "user_genome")) %>%
  separate(
    classification, 
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(across(domain:species, ~ sub(".*__", "", .)))

# ----------------------------
# Create color map for METABOLIC modules
# ----------------------------
color_map <- setNames(color_raw$color_func, color_raw$module)

color_long <- color_raw %>%
  mutate(hmm_files_singles = as.character(hmm_files), .keep = "unused") %>%
  separate_rows(hmm_files_singles, sep = "[\\s,]+")

annotation_plot <- annotation_classified %>%
  mutate(hmm_files_singles = as.character(`Hmm file`)) %>%
  separate_rows(hmm_files_singles, sep = "[\\s,]+") %>%
  left_join(color_long) %>%
  left_join(sample_key) %>%
  mutate(Category = ifelse(
    Category == "Sulfur cycling enzymes (detailed)",
    "Sulfur cycling",
    Category
  ))

# ----------------------------
# METABOLIC plot
# ----------------------------
annotation_plot %>%
  filter(!is.na(module), !is.na(isolate_id), `Hit numbers` > 1) %>%
  ggplot(aes(x = isolate_id, y = `Gene abbreviation`, color = module)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_map) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    strip.placement = "outside",
    legend.position = "bottom",
    text = element_text(size = 20),
    panel.grid.major = element_line(color = "lightgray")
  ) +
  guides(size = "none") +
  labs(x = "", y = "gene", color = "function") +
  facet_nested(rows = vars(Category), cols = vars(phylum), scales = "free", space = "free")

ggsave("results/isolate_metabolic.pdf", width = 34, height = 17, limitsize = FALSE)

# ----------------------------
# EcoFoldDB plot
# ----------------------------
eco_df %>%
  ggplot(aes(y = isolate_id, x = Gene, color = subcat)) +
  geom_point() +
  scale_color_manual(
    values = setNames(ef_color_raw$color, ef_color_raw$Subcategory)
  ) +
  labs(color = "") +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    panel.grid.major = element_line(color = "lightgray"),
    strip.placement = "outside",
    legend.position = "bottom"
  ) +
  facet_nested(cols = vars(subcat), rows = vars(phylum),
    scales = "free", space = "free") + 
  labs(x = "", y = "") + 
  guides(color = "none")

ggsave("results/ecofold_vis.png", width = 30, height = 15, limitsize = FALSE)
