## make info for corr_of_corr jobs 
library(tidyverse)
library(gmp) #### for factorize
combined.emap.data <- read_delim("basic_E-MAP_data/preprocessed_ubermap_all.txt", delim = "\t")
clusters.table <- read_delim("basic_E-MAP_data/GO_slims_2018-03-05_pearson_complete_clusters.txt", delim = "\t")
clusters <- clusters.table %>% pull(cluster) %>% unique()
#### add a randomized score column to the combined.emap.data
##### the score is ramdomized by random sampling of the scores for that Gene_uniq
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
#step <- 37*61  ### for correlations calculations
step <- 59*11*5*3   ### this is for corr of corr calculations for mut:partner pairs
pair_sets <- seq(1, n_pairs, step)
task_counter <- 1
for (i in seq_along(clusters)) {
  cluster <- clusters[i]
  for (t in seq_along(pair_sets)) {
    first_pair <- as.numeric(as.character(pair_sets[t]))
    last_pair <- first_pair + step - 1
    task_genes_and_mutants <- unique(c(pairs[1, first_pair:last_pair], pairs[2, first_pair:last_pair]))
    task.info <- list()
    task.info[["pairs"]] <- pairs[,first_pair:last_pair]
    task.info[["cluster_n"]] <- i
    task.info[["cluster"]] <- cluster
    outfilename <- str_c("clustered_correlations/info_for_mut_partner_corr_of_corr/", 
                       task_counter, "_task_info.RData" )
    save(task.info, file = outfilename)
    task_counter <- task_counter + 1
  }
}

