library(tidyverse)

# Based on this full EMAP clustering, cut into 1200 clusters, and considered
# those clusters with at least 3 members. After examining the SGD 
# descriptions, about half (21/44) strongly enrich one biological process or
# complex.

mutants_ordered_by_res <- scan('mut_orders/gsp1_mutants_ordered_by_residue.txt', what="", sep="\n")

array_descriptions_and_clusters <-
  read_delim('datasets/unbiased_array_clusters_with_descriptions.txt', delim='\t')

interesting_cluster_labels <-
  read_delim('clusters/subset_cluster_names.txt', delim=',')

interesting_clusters <-
  inner_join(array_descriptions_and_clusters,
             interesting_cluster_labels,
             by = 'cluster_number')

dir.create(file.path(getwd(), 'datasets/'))

read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t') %>% 
  gather(-mutant, key = strain, value = score) %>% 
  rename(name = strain) %>% 
  inner_join(interesting_clusters, by = 'name') %>% 
  select(-description, -name_meaning, -ORF, -cluster_number) %>% 
  group_by(cluster_name) %>% 
  nest(.key = 'cluster_data') %>% 
  mutate(cluster_data = map(cluster_data, ~spread(.x, key = name, value = score))) %>%
  mutate(cluster_data = map(cluster_data, ~mutate(.x, mutant = factor(mutant, levels = mutants_ordered_by_res)))) %>% 
  mutate(cluster_data = map(cluster_data, ~arrange(.x, mutant))) %>% 
  mutate(data = walk2(cluster_data,
                      paste0('datasets/', cluster_name, '_emap.txt'),
                      delim = '\t', na = '', write_delim))


mutants_with_GAP_GEF_data <-
  read_delim('datasets/GAP_GEF_ratio_for_clustering.txt', delim = '\t') %>% 
  pull(mutant)

read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t') %>% 
  gather(-mutant, key = strain, value = score) %>% 
  rename(name = strain) %>% 
  filter(mutant %in% mutants_with_GAP_GEF_data) %>% 
  inner_join(interesting_clusters, by = 'name') %>% 
  select(-description, -name_meaning, -ORF, -cluster_number) %>% 
  group_by(cluster_name) %>% 
  nest(.key = 'cluster_data') %>% 
  mutate(cluster_data = map(cluster_data, ~spread(.x, key = name, value = score))) %>%
  mutate(cluster_data = map(cluster_data, ~mutate(.x, mutant = factor(mutant, levels = mutants_ordered_by_res)))) %>% 
  mutate(cluster_data = map(cluster_data, ~arrange(.x, mutant))) %>% 
  mutate(data = walk2(cluster_data,
                      paste0('datasets/', cluster_name, '-GAPGEF_emap.txt'),
                      delim = '\t', na = '', write_delim))
