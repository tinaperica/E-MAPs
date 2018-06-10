library(tidyverse)
#### calculate PCA based on the correlation of correlations data
load("20180507_corr_of_corr_filter_FDR.RData")
cor_for_pca <- combined_filtered_data %>%
  mutate("geneB_cluster" = str_c(geneB, cluster, sep = " ")) %>% 
  select(geneB_cluster, geneA, corr) %>% 
  spread(geneB_cluster, corr)
cor_for_pca <- cor_for_pca[, -which(colMeans(is.na(cor_for_pca)) > 0.5)]


combined_filtered_data_pca <- cor_for_pca %>% 
  nest() %>% 
  mutate(pca = map(., ~ prcomp(.x %>% select(-geneA),
                                  center = TRUE, scale = TRUE)),
         pca_aug = map2(pca, data, ~augment(.x, data = .y)))

save(combined_filtered_data_pca, file = "20180508_pca.R")