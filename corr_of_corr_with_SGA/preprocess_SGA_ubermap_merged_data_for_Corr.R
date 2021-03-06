library(tidyverse)
library(gmp) #### for factorize

bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')

load("basic_E-MAP_data/spitzemapko.rda")
array_ORFs_from_SGA <- spitzemapko %>% 
  filter(interaction_network != "gsp1_pEMAP") %>% 
  pull(array_ORF) %>% unique()
array_orfs_to_keep <- spitzemapko %>% 
  filter(interaction_network == "gsp1_pEMAP" & ! (array_ORF %in% bad_strains) & array_ORF %in% array_ORFs_from_SGA) %>% 
  pull(array_ORF) %>% unique()

spitzemapko_for_corr <- spitzemapko %>% 
  filter(array_ORF %in% array_orfs_to_keep) %>% 
  select(query_allele_name, query_ORF, array_ORF, score)

all_queries <- spitzemapko_for_corr %>% pull(query_allele_name) %>% unique()

###### remove query_strain_ids with too many NA
queries_to_keep <- spitzemapko_for_corr %>% 
  filter(! is.na(score)) %>% 
  group_by(query_allele_name) %>% 
  summarise("count_array" = n()) %>% 
  arrange(count_array) 
queries_to_keep %>% 
  ggplot(aes(x = count_array)) + geom_density()
queries_to_keep <- queries_to_keep %>% 
  filter(count_array > 1000)
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  filter(query_allele_name %in% queries_to_keep$query_allele_name)


### add a weight to the score based on a score confidence spline
# The fit spline produces the following values:
#   -6    -5    -4    -3    -2    -1    0     1     2     3     4     5     6
# 1.00  0.99  0.95  0.83  0.43  0.14  0.02  0.07  0.26  0.65  0.86  0.95  0.99
# 
# Note that the spline fits S-score < -5.6 to values slightly (by 1e-3) above 1.
# If these confidence values are being used as weights for scores, simply cap all
# weights at 1.0.
d <- data.frame(avg_S_score = seq(-6, 6, by = 0.5),
                confidence = c(1, 1, 0.99, 0.98, 0.95, 0.9, 0.85, 0.65, 0.45,
                               0.2, 0.15, 0.08, 0, 0.05, 0.07, 0.15,
                               0.22, 0.5, 0.65, 0.75, 0.87, 0.93, .95, .98, 0.99))
confidence_spline <- smooth.spline(x = d$avg_S_score, y = d$confidence, spar = 0.3)
predict(confidence_spline, -2)
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  mutate("weight" = predict(confidence_spline, score)$y) %>% 
  mutate("weight" = ifelse(weight > 1, 1, weight))


### now make task.info files that will define which pairwise correlations are calculated in 
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

mutants <- spitzemapko_for_corr %>% 
  filter(grepl(query_allele_name, pattern = "GSP1 - ")) %>% 
  pull(query_allele_name) %>% unique()

all_genes_and_mutants <- spitzemapko_for_corr %>%
  pull(query_allele_name) %>% unique()
length(all_genes_and_mutants)
genes <- all_genes_and_mutants[! all_genes_and_mutants %in% mutants]
length(genes)

ms_hits_orfs <- read_tsv("SAINT_MS_hits.txt", col_names = F) %>%
  pull(X1)
ms_hits_gene_uniqs <- spitzemapko_for_corr %>% 
  filter(query_ORF %in% ms_hits_orfs) %>% 
  filter(! grepl(query_allele_name, pattern =  "GSP1 - ")) %>% 
  pull(query_allele_name) %>% unique()
### all the pairs, including the mutants, to calculate all correlations/similarities
all_pairs <- combn(all_genes_and_mutants, 2) ### all pairs (including mutants)
pairs_to_save <- tibble("gene1" = all_pairs[1, ], "gene2" = all_pairs[2, ]) 
write_tsv(pairs_to_save, "corr_of_corr_with_SGA/all_pairs.txt")
# only pairs of mutants and partners (for quick correlation of correlations calculations )
mut_gene_pairs <- t(as.matrix(expand.grid(mutants, genes)))
###
mut_ms_hit_pairs <- t(as.matrix(expand.grid(mutants, ms_hits_gene_uniqs)))

### make task files for all pairs
(n_pairs <- length(all_pairs)/2)
(factors <- factorize(n_pairs))
step <- 19 * 53
n_pairs / step
### 1-4564731:1007 for correlations - 4533 jobs
make_task_files(all_pairs, step, "corr_of_corr_with_SGA/all_correlations_task_info/")

#load("corr_of_corr_with_SGA/all_correlations_task_info/2129806_task_info.RData")

# # # make task file for all mut-gene pairs (for corr of corr calculations later)
# (n_pairs <- length(mut_gene_pairs)/2)
# (factors <- factorize(n_pairs))
# step <- n_pairs
# #   ## 1 job per selected cluster
# make_task_files(mut_gene_pairs, step, "corr_of_corr_with_SGA/mut_gene_corr_of_corr_task_info")
# 
# 
# # # make task file for all mut-gene pairs (for corr of corr calculations later)
# (n_pairs <- length(mut_ms_hit_pairs)/2)
# (factors <- factorize(n_pairs))
# step <- n_pairs
# #   ## 1 job per selected cluster
# make_task_files(mut_ms_hit_pairs, step, "corr_of_corr_with_SGA/mut_ms_hit_corr_of_corr_task_info/")

spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  select(-query_ORF)
#### save spitzemapko for correlations
save(spitzemapko_for_corr, file = "basic_E-MAP_data/spitzemapko_for_corr.rda")
