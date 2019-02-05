### make task_info files for pairwise correlations between Gsp1 mutants
## for all the clusters from clustered_correlations/clusters/GO_slims_2018-06-14_pearson_complete_clusters.txt"
### 3) load all the clusters
### 4) make a bunch of final .RData files (containing only info on the pairs to calculate correlations that will be calculated in that job/task) 

all_genes_and_mutants_for_pairwise_cor_calculations <- combined.emap.data %>%
  pull(Gene_uniq) %>% unique()
length(all_genes_and_mutants_for_pairwise_cor_calculations)
# real pairs
#pairs <- combn(all_genes_and_mutants_for_pairwise_cor_calculations, 2)
#### debugging pairs - just mutants
#pairs <- combn(all_genes_and_mutants_for_pairwise_cor_calculations[1:59], 2)
# only pairs of mutants and partners (useful for quick correlation of correlations calculations )
pairs <- t(as.matrix(expand.grid(
  mutants, all_genes_and_mutants_for_pairwise_cor_calculations[
    60:length(all_genes_and_mutants_for_pairwise_cor_calculations)])))

(n_pairs <- length(pairs)/2)
(factors <- factorize(length(pairs)/2))
#step <- as.numeric(as.character(factors[2]))
step <- 37*61  ### for correlations calculations
tasks <- seq(1, n_pairs, step)
## 1-10185841:2257  ### for correlations
## 1-262845:177  ## this is for mut:partner corr of corr (1485 jobs)
ubermap <- list()
ubermap[["library_clusters"]] <- clusters.table
ubermap[["clusters"]] <- clusters
ubermap[["ubermap"]] <- combined.and.randomized.emap.data
#save(ubermap, file = str_c(Sys.Date(), "_emap_data_for_corr.RData"))
task.info <- list()
for (t in seq_along(tasks)) {
  first_pair <- as.numeric(as.character(tasks[t]))
  last_pair <- first_pair + step - 1
  task_genes_and_mutants <- unique(c(pairs[1, first_pair:last_pair], pairs[2, first_pair:last_pair]))
  #task.emap.data <- combined.and.randomized.emap.data %>%
  # filter(Gene_uniq %in% task_genes_and_mutants)
  #task.info[["all_genes_and_mutants"]] <- task_genes_and_mutants
  #ubermap[["ubermap"]] <- task.emap.data
  task.info[["pairs"]] <- pairs[,first_pair:last_pair]
  outfilename <- str_c("clustered_correlations/info_for_mut_partner_corr_of_corr/", 
                       first_pair, "_task_info.RData" )
  save(task.info, file = outfilename)
}

