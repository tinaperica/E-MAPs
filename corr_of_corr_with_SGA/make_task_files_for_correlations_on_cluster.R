#### this script prepares the files for the fast and efficient calculations of correlations on the cluster
# Things it does:
### 1) defines all possible pairs of genes and mutants
### 2) make pairs and combines everything into an ubermap list used for correlation calculations on the cluster 
library(tidyverse)
library(gmp) #### for factorize
load("basic_E-MAP_data/20181230_ubermap_for_correlations.RData")

Gene_uniq_ORF_gene_name_matcher <- ubermap %>% 
  select(query_uniq, query_ORF, query_gene_name) %>% 
  unique()
write_tsv(Gene_uniq_ORF_gene_name_matcher, path = "basic_E-MAP_data/20190101_SGAubermap_gene_uniq_matcher.txt")

### now make task.info files that will define which pairwise correlations/similarities are calculated in 
## that particular task on the cluster
# a general function that makes task files when pairs and steps are defined
make_task_files <- function(pairs, step, task_outpath) {
  tasks <- seq(1, length(pairs)/2, step)
  for (t in seq_along(tasks)) {
    first_pair <- tasks[t]
    last_pair <- first_pair + step - 1
    task.info <- list()
    task.info[["pairs"]] <- pairs[, first_pair:last_pair]
    outfilename <- str_c(task_outpath, first_pair, "_task_info.RData", sep = "")
    save(task.info, file = outfilename)
  }
}

mutants <- ubermap %>% 
  filter(query_ORF == "YLR293C") %>% 
  pull(query_uniq) %>% unique()

all_genes_and_mutants <- ubermap %>%
  pull(query_uniq) %>% unique()
length(all_genes_and_mutants)
genes <- all_genes_and_mutants[! all_genes_and_mutants %in% mutants]
length(genes)

ms_hits_orfs <- read_tsv("SAINT_MS_hits.txt", col_names = F) %>%
  pull(X1)
ms_hits_gene_uniqs <- ubermap %>% 
  filter(query_ORF %in% ms_hits_orfs) %>% 
  filter(query_gene_name != "GSP1") %>% 
  pull(query_uniq) %>% unique()
### all the pairs, including the mutants, to calculate all correlations/similarities
all_pairs <- combn(all_genes_and_mutants, 2) ### all pairs (including mutants)
# only pairs of mutants and partners (for quick correlation of correlations calculations )
mut_gene_pairs <- t(as.matrix(expand.grid(mutants, genes)))
###
mut_ms_hit_pairs <- t(as.matrix(expand.grid(mutants, ms_hits_gene_uniqs)))

### make task files for all pairs
(n_pairs <- length(all_pairs)/2)
(factors <- factorize(n_pairs))
step <- 41 * 53
### 1-21245421:2173 for correlations - 9777 jobs
make_task_files(all_pairs, step, "corr_of_corr_with_SGA/all_correlations_task_info/")

# # make task file for all mut-gene pairs (for corr of corr calculations later)
(n_pairs <- length(mut_gene_pairs)/2)
(factors <- factorize(n_pairs))
step <- n_pairs
#   ## 1 job per selected cluster
make_task_files(mut_gene_pairs, step, "corr_of_corr_with_SGA/mut_gene_corr_of_corr_task_info")


# # make task file for all mut-gene pairs (for corr of corr calculations later)
(n_pairs <- length(mut_ms_hit_pairs)/2)
(factors <- factorize(n_pairs))
step <- n_pairs
#   ## 1 job per selected cluster
make_task_files(mut_ms_hit_pairs, step, "corr_of_corr_with_SGA/mut_ms_hit_corr_of_corr_task_info/")
load("corr_of_corr_with_SGA/mut_ms_hit_corr_of_corr_task_info/1_task_info.RData")

ubermap <- data %>% 
  select(Gene_uniq, library.ORF, score)
save(ubermap, file = "basic_E-MAP_data/20180921_emap_data_for_corr.RData")

