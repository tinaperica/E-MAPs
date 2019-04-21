library(tidyverse)
set.seed(94158)

# Cut EMAP into k clusters, and consider only those clusters with at least 3
# members. After examining the SGD descriptions, about half (21/44) strongly
# enrich one biological process or complex. These clusters were manually labeled
# and the file 'subset_cluster_names.txt' was made.
# NOTE: different clusters if you filter to GAP/GEF complete first

# read in dataset, as matrix, then cluster
mat <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
dissim <- as.dist((1 - cormat)/2)
hc <- hclust(dissim, method = "average" )

k = 1200 # set k = number of clusters to make
cluster_assignments <- cutree(hc, k = k)

clusters_with_descriptions <-
  data.frame('name' = names(cluster_assignments), 'cluster_number' = unname(cluster_assignments), stringsAsFactors = FALSE) %>%
  left_join(read_delim('SGD_descriptions.txt', delim='\t', col_types = cols()), by = 'name')

write_delim(clusters_with_descriptions, 'unbiased_array_clusters_with_descriptions.txt', delim='\t')

cluster_dir = 'clusters/'
dir.create(file.path(getwd(), cluster_dir))

clusters_with_descriptions %>%
  group_by(cluster_number) %>%
  filter(n() > 2) %>%
  nest(.key = 'cluster_data') %>%
  mutate(data = walk2(cluster_data,
                      paste0(cluster_dir, '/', cluster_number, '.csv'),
                      delim = ',', na = '', write_delim))

gene_sets_from_clustering <-
  clusters_with_descriptions %>%
  inner_join(read_delim('subset_cluster_names.txt', delim = ',', col_types = cols())) %>%
  arrange(cluster_number) %>%
  select(ORF, name, cluster_name) %>%
  rename('Gene' = name, 'gene_set'  = cluster_name) %>%
  mutate('set_label' = 'from_clustering')

gene_sets_processes <-
  read_delim('gene_sets_for_gsp1_functions.txt', delim='\t', col_types = cols()) %>%
  rename('gene_set' = term) %>%
  mutate('set_label' = 'curated_processes')

gene_sets_complexes <-
  read_delim('all_wodak_complexes.txt', delim='\t', col_types = cols()) %>%
  rename('Gene' = Name, 'gene_set' = Complex) %>%
  mutate('set_label' = 'wodak_complexes')

gene_sets <- bind_rows(gene_sets_from_clustering, gene_sets_processes, gene_sets_complexes)

write_delim(gene_sets, 'gene_sets_combined.txt', delim='\t')

gene_sets_min3 <-
  bind_rows(gene_sets_from_clustering, gene_sets_processes, gene_sets_complexes) %>%
  group_by(gene_set) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count > 2) %>%
  select(-count)

write_delim(gene_sets_min3, 'gene_sets_combined_min3.txt', delim='\t')