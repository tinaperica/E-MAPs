#### this script prepares the files for the fast and efficient calculations of correlations on the cluster
# Things it does:
### 1) defines all possible pairs of genes and mutants
### 2) makes random E-MAP scores - for each query: same scores but library genes names will be shuffled
### 3) randomizes only the high scoring library genes
### 4) make pairs and combine everything into an ubermap list used for correlation calculations on the cluster 
library(tidyverse)
library(gmp) #### for factorize
data <- read_tsv("basic_E-MAP_data/preprocessed_ubermap_all.txt") # includes data for Gsp1 mutants and the ubermap data
clusters.file <- "clustered_correlations/clusters/GO_slims_2018-06-30_pearson_complete_clusters.txt"
clusters.table <- read_tsv(clusters.file)   # all the library genes GO and GI based categories
selected_clusters_file <- "choose_clusters/2018-07-02_selected_pearson_complete_clusters.txt"
if ( file.exists(selected_clusters_file)) {
  selected.clusters.table <- read_tsv(selected_clusters_file) 
} else {
  selected.clusters.table <- read_tsv(clusters.file)
}
high_score_lib_genes <- read_tsv("library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt")
clusters <- clusters.table %>% pull(cluster) %>% unique()
selected.clusters <- selected.clusters.table %>% pull(cluster) %>% unique()
#### on 20180917 - added extended selected clusters
extended.selected.clusters <- sort(clusters[grepl(clusters, pattern = "GO_30_[0-9]+_mut", perl = T)])
extendend.selected.clusters <- sort(unique(append(extended.selected.clusters, selected.clusters)))
#### add a randomized score column to the combined.emap.data
##### the score is ramdomized by random sampling of the scores for that Gene_uniq
### two types of randomization: all and high
set.seed(2013)
## first shuffle only the scores of the high scoring library genes (defined in the high_score_lib_genes)
high_score_randomized_ubermap <- data %>% 
  filter(library.ORF %in% high_score_lib_genes$library.ORF) %>% 
  group_by(Gene_uniq) %>% 
  mutate("random_high_score" = sample(x = score, size = length(score)))
randomized_ubermap <- data %>% 
  filter(! library.ORF %in% high_score_randomized_ubermap$library.ORF) %>% 
  mutate("random_high_score" = score) %>% 
  bind_rows(., high_score_randomized_ubermap) %>% 
  arrange(Gene_uniq, library.ORF)
### then shuffle all the scores
randomized_ubermap <- randomized_ubermap %>%
  group_by(Gene_uniq) %>%
  mutate("random_score" = sample(x = score, size = length(score)))

### now do everything once more with a different seed
set.seed(2017)
high_score_randomized_ubermap <- randomized_ubermap %>% 
  filter(library.ORF %in% high_score_lib_genes$library.ORF) %>% 
  group_by(Gene_uniq) %>% 
  mutate("random_high_score_2" = sample(x = score, size = length(score)))
randomized_ubermap <- randomized_ubermap %>% 
  filter(! library.ORF %in% high_score_randomized_ubermap$library.ORF) %>% 
  mutate("random_high_score_2" = score) %>% 
  bind_rows(., high_score_randomized_ubermap) %>% 
  arrange(Gene_uniq, library.ORF)
### shuffle all the scores again
randomized_ubermap <- randomized_ubermap %>%
  group_by(Gene_uniq) %>%
  mutate("random_score_2" = sample(x = score, size = length(score)))

### now make task.info files that will define which pairwise correlations/similarities are calculated in 
## that particular task on the cluster
# a general funciton that makes task files when pairs and steps are defined
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

mutants <- data %>% 
  filter(ORF == "YLR293C") %>% 
  pull(Gene_uniq) %>% unique()

all_genes_and_mutants <- data %>%
  pull(Gene_uniq) %>% unique()
length(all_genes_and_mutants)
genes <- all_genes_and_mutants[! all_genes_and_mutants %in% mutants]
length(genes)

### only mutant pairs (for choosing which library (sub)clusters are informative)
cluster_eval_pairs <- combn(mutants, 2) ### only mutant pairs
### all the pairs, including the mutants, to calculate all correlations/similarities
all_pairs <- combn(all_genes_and_mutants, 2) ### all pairs (including mutants)
# only pairs of mutants and partners (for quick correlation of correlations calculations )
mut_gene_pairs <- t(as.matrix(expand.grid(mutants, genes)))

### first make task.info files for mutant only pairwise calculations
(n_pairs <- length(cluster_eval_pairs)/2)
(factors <- factorize(n_pairs))  
step <- as.numeric(as.character(factors[1]))
#### 1-1711:29
make_task_files(pairs = cluster_eval_pairs, step, "clustered_correlations/mutant_correlations_task_info/")

### now for all pairs
(n_pairs <- length(all_pairs)/2)
(factors <- factorize(n_pairs))
step <- 181  * 5
### 1-10240075:905  ## for correlations - 11315 jobs
make_task_files(all_pairs, step, "clustered_correlations/all_correlations_task_info/")

# # finally for all mut-gene pairs (for corr of corr calculations later)
(n_pairs <- length(mut_gene_pairs)/2)
(factors <- factorize(n_pairs))
step <- n_pairs
# 1-183  ## 1 job per selected cluster
make_task_files(mut_gene_pairs, step, "clustered_correlations/mut_gene_corr_of_corr_task_info/")

ubermap <- list()
ubermap[["ubermap"]] <- randomized_ubermap
ubermap[["library_clusters"]] <- clusters.table
ubermap[["clusters"]] <- clusters
save(ubermap, file = "basic_E-MAP_data/20180630_emap_data_for_corr_all_clusters.RData")

ubermap <- list()
ubermap[["ubermap"]] <- randomized_ubermap
ubermap[["library_clusters"]] <- selected.clusters.table
ubermap[["clusters"]] <- selected.clusters
save(ubermap, file = "basic_E-MAP_data/20180630_emap_data_for_corr.RData")

