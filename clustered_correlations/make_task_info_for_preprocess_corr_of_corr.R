library(tidyverse)
pair_sets <- 1:27
#pair_sets <- 1:2
smaller_pair_sets <- seq(1, 9735, 177)
file_counter <- 1
corr_of_corr_file_path <- "~/E-MAP_analysis/Output/corr_of_corr_mut_per_cluster/20180409_correlations_of_correlations_mut_partner_pairs_only_"
#corr_of_corr_file_path <- "clustered_correlations/temp/20180409_correlations_of_correlations_mut_partner_pairs_only_"

for (ps in pair_sets) {
  file_numbers <- seq(ps, (27*279), 27)
  all_files_per_task_set <- tibble(geneA = character(), geneB = character(), 
                                   fn = integer(), cluster = character(), cluster_n = integer())
  for (fn in file_numbers) {
    file <- str_c("clustered_correlations/info_for_mut_partner_corr_of_corr/", fn, "_task_info.RData", collapse = "")
    load(file)
    cluster_n <- task.info[["cluster_n"]]
    cluster <- task.info[["cluster"]]
    geneA <- task.info[["pairs"]][1, 1:length(task.info[["pairs"]][1,])]
    geneB <- task.info[["pairs"]][2, 1:length(task.info[["pairs"]][1,])]
    all_files_per_task_set <-  bind_rows(all_files_per_task_set,
                    tibble(geneA, geneB, fn = as.integer(fn), cluster, cluster_n))
  }
  all_files_per_task_set <- all_files_per_task_set %>% 
    mutate(pair = str_c(geneA, geneB, sep = " "), 
           filename = str_c(corr_of_corr_file_path, cluster_n, "_", fn, ".RData" )) 
  rm(task.info)
  pairs <- all_files_per_task_set %>% pull(pair) %>% unique()
  n_pairs <- length(pairs) ### should be always 9735
  for (sps in smaller_pair_sets) {
    first_pair <- sps
    last_pair <- first_pair + 177 - 1
    temp_pairs <- pairs[first_pair:last_pair]
    task.info <- all_files_per_task_set %>% 
      filter(pair %in% temp_pairs) %>%
      select(-pair)
      
    outfile <- str_c("clustered_correlations/info_for_corr_of_corr_processing/", file_counter, "_task_info.RData")
    save(task.info, file = outfile)
    file_counter <- file_counter + 1
  }
}

