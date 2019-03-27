library(tidyverse)
load("basic_E-MAP_data/spitzemapko.rda")
query_orf_index <- spitzemapko %>% 
  select(query_allele_name, query_ORF) %>% unique()
#pairs <- read_tsv("corr_of_corr_with_SGA/all_pairs.txt")
corr_path <- "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/March2019_analysis/correlations/"
(files <- dir(corr_path))
files <- files[grepl(files, pattern = ".RData")]

for (f in seq_along(files)) {
    file_path <- file.path(corr_path, files[f])
    print(file_path)
    load(file_path)
    correlations_df <- as_tibble(correlations_df)
    outfilepath <- file.path("corr_of_corr_with_SGA/correlations", str_c(files[f], ".txt"))
    write_tsv(correlations_df, path = outfilepath)
}

corr_path <- "corr_of_corr_with_SGA/correlations/"
(files <- dir(corr_path))
files <- files[grepl(files, pattern = ".RData.txt")]

correlations <- files %>%
  map(function(x) read_tsv(file.path(corr_path, x))) %>%  
  reduce(bind_rows)

output_file <- str_c("corr_of_corr_with_SGA/spitzemapko_correlations_all.RData")
save(correlations, file = output_file)

# pearson_correlations <- correlations %>% 
#   select(query_uniq1, query_uniq2, pearson)
# min_w_pearson <- correlations %>% 
#   select(query_uniq1, query_uniq2, min_weighted_pearson)
# max_w_pearson <- correlations %>% 
#   select(query_uniq1, query_uniq2, max_weighted_pearson)
# mean_w_pearson <- correlations %>% 
#   select(query_uniq1, query_uniq2, mean_weighted_pearson)
# 
# 
# 
# output_file <- str_c("corr_of_corr_with_SGA/correlations/spitzemapko_correlations_pearson.RData")
# save(pearson_correlations, file = output_file)
# 
# output_file <- str_c("corr_of_corr_with_SGA/correlations/spitzemapko_correlations_min_weighted.RData")
# save(min_w_pearson, file = output_file)
# 
# output_file <- str_c("corr_of_corr_with_SGA/correlations/spitzemapko_correlations_max_weighted.RData")
# save(max_w_pearson, file = output_file)
# 
# output_file <- str_c("corr_of_corr_with_SGA/correlations/spitzemapko_correlations_mean_weighted.RData")
# save(mean_w_pearson, file = output_file)

Gsp1_mutants_correlations <- correlations %>% 
  filter(grepl(query_uniq1, pattern = "GSP1 - ") & ! grepl(query_uniq2, pattern = "GSP1 - "))
output_file <- str_c("corr_of_corr_with_SGA/correlations/spitzemapko_all_correlations_mutants_only.RData")
save(Gsp1_mutants_correlations, file = output_file)


### now format it for gia
### make left annotation
mutants <- Gsp1_mutants_correlations %>% 
  select(query_uniq1) %>% 
  unique() %>% 
  separate(query_uniq1, into = c("gsp1", "mut"), remove = F, sep = " - ") %>% 
  select(-gsp1)
write_tsv(mutants, path = "gia_analysis/corr_Gsp1_mut_left_annotation.txt", col_names = F)
### make the right annotation for gia (gene sets with unique SGA names)
gene_sets <- read_tsv("GSEA_like_analysis/gene_groups.txt")
right_annotations_correlations <- Gsp1_mutants_correlations %>% 
  select(query_uniq2) %>% 
  unique() %>% 
  inner_join(., query_orf_index, by = c("query_uniq2" = "query_allele_name")) %>% 
  inner_join(., gene_sets, by = c("query_ORF" = "ORF")) %>% 
  select(query_uniq2, term)
write_tsv(right_annotations_correlations, path = "gia_analysis/corr_Gsp1_mut_right_annotation.txt", col_names = F)

#### make the spread data for gia
Gsp1_mutants_correlations <- Gsp1_mutants_correlations %>% 
  select(query_uniq1, query_uniq2, pearson) %>% 
  spread(query_uniq2, pearson)
write_tsv(Gsp1_mutants_correlations, "gia_analysis/Gsp1_spitzemapko_correlations_for_gia.txt", na = "")



#### output correlations for correlations of correlations
## that means remove all cases of query_allele_name_2 being a GSP1 mutant from our screen
correlations <- correlations %>% 
  filter(! grepl(query_uniq2, pattern = "GSP1 - "))

corr <- correlations %>% 
  select(query_uniq1, query_uniq2, 'corr' = pearson)
save(corr, file = "corr_of_corr_with_SGA/correlations_for_corr_of_corr/pearson_for_corr_of_corr.RData")

corr <- correlations %>% 
  select(query_uniq1, query_uniq2, 'corr' = min_weighted_pearson)
save(corr, file = "corr_of_corr_with_SGA/correlations_for_corr_of_corr/min_pearson_for_corr_of_corr.RData")

corr <- correlations %>% 
  select(query_uniq1, query_uniq2, 'corr' = max_weighted_pearson)
save(corr, file = "corr_of_corr_with_SGA/correlations_for_corr_of_corr/max_pearson_for_corr_of_corr.RData")

corr <- correlations %>% 
  select(query_uniq1, query_uniq2, 'corr' = mean_weighted_pearson)
save(corr, file = "corr_of_corr_with_SGA/correlations_for_corr_of_corr/mean_pearson_for_corr_of_corr.RData")

