#### this script prepares the files for the fast and efficient calculations of correlations on the cluster
# Things it does:
### 1) defines all possible pairs of genes and mutants
### 2) makes random E-MAP scores - for each query: same scores but library genes names will be shuffled
### 3) load all the clusters
### 4) make a bunch of final .RData files (containing only info on the pairs to calculate correlstions that will be calulted in that job/task) 
library(tidyverse)
library(gmp) #### for factorize
combined.emap.data <- read_delim("basic_E-MAP_data/preprocessed_ubermap_all.txt", delim = "\t")
clusters.table <- read_delim("basic_E-MAP_data/GO_slims_2018-03-05_pearson_complete_clusters.txt", delim = "\t")
clusters <- clusters.table %>% pull(cluster) %>% unique()
#### add a randomized score column to the combined.emap.data
##### the score is ramdomized by random sampling of the scores for that Gene_uniq
set.seed(2013)
combined.and.randomized.emap.data <- combined.emap.data %>%
  group_by(Gene_uniq) %>%
  mutate("random_score" = sample(x = score, size = length(score)))
set.seed(2017)
combined.and.randomized.emap.data <- combined.and.randomized.emap.data %>%
  group_by(Gene_uniq) %>%
  mutate("random_score_2" = sample(x = score, size = length(score)))
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","NTER3XFLAG WT","CTER3XFLAG WT","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T")
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

