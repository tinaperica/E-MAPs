#### this script prepares the files for the fast and efficient calculations of correlations on the cluster
# Things it does:
### 1) defines all possible pairs of genes and mutants
### 2) makes a random version of a preprocessed_ubermap_all.txt which will have the same scores but library genes names will be shuffled
### 3) load all the clusters
### 4) make a bunch of final .RData files (containing only data for the files that will be calulted in that job/task) 
library(tidyverse)
library(gmp) #### for factorize
combined.emap.data <- read_delim("basic_E-MAP_data/preprocessed_ubermap_all.txt", delim = "\t")
clusters.table <- read_delim("basic_E-MAP_data/GO_slims_2018-03-05_pearson_complete_clusters.txt", delim = "\t")
clusters <- clusters.table %>% pull(cluster) %>% unique()
#### add a randomized score column to the combined.emap.data
##### the score is ramdomized by random sampling of the scores for that Gene_uniq
combined.and.randomized.emap.data <- combined.emap.data %>%
  group_by(Gene_uniq) %>%
  mutate("random_score" = sample(x = score, size = length(score)))
##### test plots
# combined.and.randomized.emap.data %>% 
#   filter(Gene_uniq == "F58A") %>%
#   ggplot(aes(x = score)) + geom_density() + geom_density(aes(x = random_score), color = "red", linetype = 2)
# combined.and.randomized.emap.data %>% 
#   filter(Gene_uniq == "R108L") %>%
#   ggplot(aes(x = score)) + geom_density() + geom_density(aes(x = random_score), color = "red", linetype = 2)

########## make all possible pairs of Gsp1 mutants and ubermap genes for correlation calculations
all_genes_and_mutants_for_pairwise_cor_calculations <- combined.emap.data %>%
  pull(Gene_uniq) %>% unique()
length(all_genes_and_mutants_for_pairwise_cor_calculations)
# real pairs
pairs <- combn(all_genes_and_mutants_for_pairwise_cor_calculations, 2)
#### debugging pairs
#pairs <- combn(all_genes_and_mutants_for_pairwise_cor_calculations[1:59], 2)

(n_pairs <- length(pairs)/2)
(factors <- factorize(length(pairs)/2))
step <- as.numeric(as.character(factors[3]))
#step <- 1711
tasks <- seq(1, n_pairs, step)

ubermap <- list()
ubermap[["library_clusters"]] <- clusters.table
ubermap[["clusters"]] <- clusters
for (t in seq_along(tasks)) {
  first_pair <- as.numeric(as.character(tasks[t]))
  last_pair <- first_pair + step - 1
  task_genes_and_mutants <- unique(c(pairs[1, first_pair:last_pair], pairs[2, first_pair:last_pair]))
  task.emap.data <- combined.and.randomized.emap.data %>%
    filter(Gene_uniq %in% task_genes_and_mutants)
  ubermap[["all_genes_and_mutants"]] <- task_genes_and_mutants
  ubermap[["ubermap"]] <- task.emap.data
  ubermap[["pairs"]] <- pairs[,first_pair:last_pair]
  outfilename <- str_c("preprocessed_data_for_correlations/", 
                       first_pair, "_emap_data_for_corr.RData" )
  save(ubermap, file = outfilename)
}
