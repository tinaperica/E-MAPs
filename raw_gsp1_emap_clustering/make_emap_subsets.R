
# Based on this full EMAP clustering, cut into 1200 clusters, and considered
# those clusters with at least 3 members. After examining the SGD 
# descriptions, about half (21/44) strongly enrich one biological process or
# complex.

library(tidyverse)

mutants_ordered_by_res <- scan('gsp1_mutants_ordered_by_residue.txt', what="", sep="\n")

array_descriptions_and_clusters <-
  read_delim('unbiased_array_clusters_with_descriptions.txt', delim='\t')

interesting_cluster_labels <-
  read_delim('subset_cluster_names.txt', delim=',')

interesting_clusters <-
  inner_join(array_descriptions_and_clusters,
             interesting_cluster_labels,
             by = 'cluster')

cluster_emaps_dir = 'clustered_emaps/'
dir.create(file.path(getwd(), cluster_emaps_dir))

read_delim('gsp1_emap_for_clustering.txt', delim = '\t') %>% 
  gather(-mutant, key = strain, value = score) %>% 
  rename(name = strain) %>% 
  inner_join(interesting_clusters, by = 'name') %>% 
  select(-description, -name_meaning, -ORF, -cluster) %>% 
  group_by(cluster_name) %>% 
  nest(.key = 'cluster_data') %>% 
  mutate(cluster_data = map(cluster_data, ~spread(.x, key = name, value = score))) %>%
  mutate(cluster_data = map(cluster_data, ~mutate(.x, mutant = factor(mutant, levels = mutants_ordered_by_res)))) %>% 
  mutate(cluster_data = map(cluster_data, ~arrange(.x, mutant))) %>% 
  mutate(data = walk2(cluster_data,
                      paste0(cluster_emaps_dir, '/', cluster_name, '_emap.txt'),
                      delim = '\t', na = '', write_delim))
