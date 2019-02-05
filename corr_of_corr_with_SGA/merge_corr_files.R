library(tidyverse)
matcher <- read_tsv("basic_E-MAP_data/20190101_SGAubermap_gene_uniq_matcher.txt")
gsp1_mutants <- matcher %>% filter(query_gene_name == "GSP1") %>% 
  pull(query_uniq) %>% unique()
corr_path <- "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Jan2019_analysis/correlations/"
(files <- dir(corr_path))
files <- files[grepl(files, pattern = ".RData")]
input_counter <- 0
output_counter <- 0
merged_corr <- data.frame()
for (f in seq_along(files)) {
  file_path <- file.path(corr_path, files[f])
  load(file_path)
  correlations_df <- as_tibble(correlations_df)
  merged_corr <- rbind(merged_corr, correlations_df)
  correlations_df <- correlations_df %>% 
    select("query_uniq1" = query_uniq2, "query_uniq2" = query_uniq1, w_correlation, soft_cos_sim)
  merged_corr <- rbind(merged_corr, correlations_df)
  input_counter <- input_counter + 1
  if (input_counter == 977) {
    output_counter <- output_counter + 1
    output_file <- str_c("corr_of_corr_with_SGA/correlations/correlations_", output_counter, ".RData")
    merged_corr <- merged_corr %>% 
      filter(! query_uniq2 %in% gsp1_mutants)
    save(merged_corr, file = output_file)
    rm(merged_corr)
    merged_corr <- data.frame()
    input_counter <- 0
  }
}
output_file <- str_c("corr_of_corr_with_SGA/correlations/correlations_10.RData")
save(merged_corr, file = output_file)
rm(merged_corr)

SGA_correlations <- tibble(query_uniq1 = character(), query_uniq2 = character(), 
                           w_correlation = double(), soft_cos_sim = double())

corr_path <- "corr_of_corr_with_SGA/correlations/"
(files <- dir(corr_path))
files <- files[grepl(files, pattern = ".RData")]
for (f in seq_along(files)) {
  file_path <- file.path(corr_path, files[f])
  load(file_path)
  SGA_correlations <- SGA_correlations %>% bind_rows(., merged_corr)
}
save(SGA_correlations, file = "corr_of_corr_with_SGA/20190102_all_correlations.RData")
