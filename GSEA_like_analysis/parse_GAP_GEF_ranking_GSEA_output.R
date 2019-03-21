#### for parsing and analysing the output from GSEA 3.0
library(tidyverse)
input <- "20190311"
input_dir <- file.path("GSEA_like_analysis/SGSEA 3.0 output folder", input)

gsea_out_folders <- dir(input_dir)
files_tib <- tibble("rank" = as.character(), "filepaths" = as.character(), "enrichment_sign" = as.character())
for (o in seq_along(gsea_out_folders)) {
  rank <- gsea_out_folders[o] %>% str_split("_") %>% .[[1]]  %>% .[1]
  dir_path <- file.path(input_dir, gsea_out_folders[o])
  split_dir <- str_split(gsea_out_folders[o], pattern = ".GseaPreranked.")
  filenames <- dir(dir_path)
  filenames <- filenames[grepl(filenames, pattern = "gsea_report_for_na.+xls", perl = T)]
  enrichment_sign <- filenames %>% str_split("_") %>% unlist() %>% .[c(5, 11)]
  filepaths <- file.path(dir_path, filenames)
  files_tib <- bind_rows(files_tib, tibble("rank" = rank, "filepaths" = filepaths, 
                        "enrichment_sign" = enrichment_sign))
}
enrichment_data <- tibble("gene_set" = as.character(), "gene_set_size"= integer(), "norm_ES" = double(), 
                          "FDR" = double(), "rank" = character(), "sign" = character())
for (i in seq_along(files_tib$filepaths))  {
  filepath <- files_tib$filepaths[i]
  rank <- files_tib$rank[i]
  enrichment_sign <- files_tib$enrichment_sign[i]
  input_tib <- read_tsv(filepath, col_names = T) %>% 
    select("gene_set" = NAME, "gene_set_size" = SIZE, 'norm_ES' = NES, "FDR" = `FDR q-val`) %>% 
    mutate("rank" = rank, "sign" = enrichment_sign)
  enrichment_data <- bind_rows(enrichment_data, input_tib) %>% 
    arrange(gene_set)
}


enrichment_diff <- enrichment_data %>% 
  group_by(gene_set) %>% 
  mutate("NES_diff" = max(norm_ES) - min(norm_ES)) %>% 
  mutate("direction" = ifelse(NES_diff < abs(norm_ES), "same", "opposite")) %>% 
  #filter(direction == "opposite") %>% 
  mutate("GAP or GEF enriched" = ifelse((sign == "pos" & rank == "GAP"), "GAP", "GEF")) %>% 
  ungroup()
gene_sets_by_enrich_diff <- enrichment_diff %>% arrange(desc(NES_diff)) %>% 
  pull(gene_set) %>% unique()
enrichment_diff %>% 
  filter(sign == "pos") %>% 
  select(gene_set, NES_diff, `GAP or GEF enriched`) %>% unique() %>% 
  filter(NES_diff >= 2.5) %>% 
  mutate("gene_set" = factor(gene_set, gene_sets_by_enrich_diff)) %>% 
  ggplot(aes(x = gene_set, y = NES_diff, fill = `GAP or GEF enriched`)) + 
  geom_bar(stat = "identity")  +
  ylab("Difference in NES (size normalized GSEA enrichment score)") + xlab("Gene set") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = margin(1, 1, 1, 3.5, "cm"))
ggsave("GSEA_like_analysis/NES_diff.pdf", height = 15, width = 25)  

enrichment_diff %>% 
  filter(`GAP or GEF enriched` == "GAP") %>% 
  select(gene_set, NES_diff) %>% unique() %>% 
  filter(NES_diff > 2) %>% 
  mutate("gene_set" = factor(gene_set, gene_sets_by_enrich_diff)) %>% 
  ggplot(aes(x = gene_set, y = NES_diff)) + 
  geom_bar(stat = "identity")  +
  ylab("Difference in NES (size normalized GSEA enrichment score)") + xlab("Gene set") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(1, 1, 1, 3.9, "cm"))
ggsave("GSEA_like_analysis/GAP_enriched_NES_diff.pdf", height = 13, width = 20)  

enrichment_diff %>% 
  filter(`GAP or GEF enriched` == "GEF") %>% 
  select(gene_set, NES_diff) %>% unique() %>% 
  filter(NES_diff > 2) %>% 
  mutate("gene_set" = factor(gene_set, gene_sets_by_enrich_diff)) %>% 
  ggplot(aes(x = gene_set, y = NES_diff)) + 
  geom_bar(stat = "identity")  +
  ylab("Difference in NES (size normalized GSEA enrichment score)") + xlab("Gene set") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(1, 1, 1, 4.2, "cm")
        )
ggsave("GSEA_like_analysis/GEF_enriched_NES_diff.pdf", height = 13, width = 23)  

gene_sets_with_good_FDR <- enrichment_diff %>% 
  filter(FDR < 0.05) %>% 
  mutate("gene_set" = ifelse(gene_set == "SAC AND APC", "SPINDLE ASSEMBLY CHECKPOINT AND APC", gene_set)) %>% 
  pull(gene_set) %>% unique()
GAP_nes_order <- enrichment_diff %>% 
  mutate("gene_set" = ifelse(gene_set == "SAC AND APC", "SPINDLE ASSEMBLY CHECKPOINT AND APC", gene_set)) %>% 
  filter(rank == "GAP" & gene_set %in% gene_sets_with_good_FDR & direction == "opposite") %>% 
  arrange(norm_ES) %>% pull(gene_set)
enrichment_diff %>% 
  mutate("gene_set" = ifelse(gene_set == "SAC AND APC", "SPINDLE ASSEMBLY CHECKPOINT AND APC", gene_set)) %>% 
  filter(gene_set %in% GAP_nes_order & direction == "opposite") %>%   
  mutate("gene_set" = factor(gene_set, GAP_nes_order)) %>% 
  ggplot(aes(x = rank, y = gene_set, fill = norm_ES, size = FDR)) + 
  scale_size("FDR", range = c(6, 0.5), breaks = c(0, 0.0375, 0.05, 0.5)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() +
  ylab("Gene set") + xlab("Enrichment in rank correlation with GAP or GEF") +
  labs(fill = "normalized\nenrichment score")
ggsave("GSEA_like_analysis/NES_point_plot.pdf", height = 3.5, width = 8)



