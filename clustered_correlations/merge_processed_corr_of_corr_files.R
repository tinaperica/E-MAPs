### merge all the processed corr_of_corr files
library(tidyverse)
library(broom)
library(ggdendro)
orf_to_protein_name_index <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F) %>% 
  select("ORF" = X1, "gene_name" = X2) %>% unique()
combined_data <- tibble()
file_names <- dir("clustered_correlations/processed_corr_of_corr/")
#file_names <- file_names[1:100]
for (fn in seq_along(file_names)) {
  file_path = file.path("clustered_correlations/processed_corr_of_corr/", file_names[fn])
  load(file_path)
  all_corr_of_corr <- as_tibble(all_corr_of_corr) %>% 
    select(cluster, geneA, geneB, corr, FDR, random_avrg_corr, random_avrg_FDR, corr_ratio) %>% 
    mutate("ORF" = gsub(x = geneB, pattern = "_.+", perl = T, replacement = "")) %>% 
    inner_join(., orf_to_protein_name_index, by = "ORF")
  combined_data <- bind_rows(combined_data, all_corr_of_corr)
}
rm(all_corr_of_corr)

save(combined_data, file = "20180507_corr_of_corr_all.RData")
 ### filter data by FDR
combined_filtered_data <- combined_data %>% 
  filter(FDR < 0.05 & random_avrg_FDR > 0.05)
save(combined_filtered_data, file = "20180507_corr_of_corr_filter_FDR.RData")
 
