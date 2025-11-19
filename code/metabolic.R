library("tidyverse")
library("ggh4x")

# set file paths
annotation_fp <- "./data/metabolic_g/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv"
color_fp <- "/projects/p32449/color_palettes/color_dictionary_draft07.csv"
classification_fp <- "./data/gtdbtk_all/gtdbtk.bac120.summary.tsv"
sample_key_fp <- "./data/sharepoint_data/submission_key.csv"

# read in data
annotation_raw <- read_tsv(annotation_fp)
color_raw <- read_csv(color_fp)
classification_raw <- read_tsv(classification_fp)
sample_key <- read_csv(sample_key_fp)

# tidy data for ggplot visualization
annotation_tidy <- annotation_raw %>% 
  mutate(across(-c(Category:`Hmm detecting threshold`), as.character)) %>% 
  pivot_longer(
    cols = -c(Category:`Hmm detecting threshold`),
    names_to = "bin_measure",
    values_to = "value") %>% 
  separate(bin_measure, into = c("sample_bin", "measurement"), 
    sep = " ",
    extra = "merge") %>% 
  pivot_wider(names_from = "measurement",
    values_from = "value") 

# add taxonomy
annotation_classified <- annotation_tidy %>% 
  mutate(sample = str_remove_all(sample_bin, "_scaffolds")) %>% 
  left_join(classification_raw, by = c("sample" = "user_genome")) %>% 
  separate(classification, into = 
    c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";", fill = "right", remove = FALSE) %>%
  mutate(across(domain:species, ~ sub(".*__", "", .)))

# create color map for plotting
color_map <- 
  setNames(color_raw$color_func, color_raw$module)

# add colors to annotation file
color_long <- color_raw %>% 
  mutate(hmm_files_singles = as.character(hmm_files), 
    .keep = "unused") %>% 
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

# visualize!
annotation_plot %>% 
  filter(!is.na(module)) %>% 
  filter(`Hit numbers` > 0, !is.na(isolate_id)) %>% 
  ggplot(aes(x = isolate_id, y = `Gene abbreviation`, color = module)) + 
    geom_point(size = 3) + 
    scale_color_manual(values = color_map) + 
    theme_gray() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          strip.text.y = element_text(angle = 0),
          strip.text.x = element_text(angle = 90), 
          strip.placement = "outside",
          legend.position = "bottom",
          text = element_text(size = 20)) + 
    guides(size = "none") + 
    labs(x = "", 
         y = "gene",
         color = "function") + 
    facet_nested(rows = vars(Category), 
                 cols = vars(phylum),
                 scales = "free", space = "free")

ggsave("results/isolate_metabolic.pdf",
  width = 34,
  height = 17, 
  limitsize = FALSE)