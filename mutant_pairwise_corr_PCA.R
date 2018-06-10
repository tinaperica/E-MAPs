library(tidyverse)
library(mclust)
library(ggcorrplot)
library(pcaMethods)

##### correlations between mutants
load("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/clustered_correlations/correlation_RData/mutants_only.RData")
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","NTER3XFLAG WT","CTER3XFLAG WT","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T")
selected_mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I")
cluster_annotation <- read_tsv("clustered_correlations/clusters/GO_slims_2018-03-05_pearson_complete_clusters.txt", col_names = T)
correlations_df <- as_tibble(correlations_df)
clusters_to_consider <- read_tsv("choose_clusters/clusters_to_consider.txt")
representative_clusters <- read_tsv("choose_clusters/representative_clusters.txt")
correlations_df_to_consider <- correlations_df %>% 
  filter(cluster %in% clusters_to_consider$clusters_to_consider)

# correlation between clusters
corr_spread <- correlations_df_to_consider %>% 
  mutate("pair" = str_c(Gene_uniq1, Gene_uniq2, sep = " ")) %>% 
  select(pair, cluster, correlation) %>%
  spread(cluster, correlation)
corr_spread_matrix <- as.matrix(corr_spread[, -1])
rownames(corr_spread_matrix) <- corr_spread$pair
##### correlations between mutants
cormat <- cor(corr_spread_matrix[,], use = "pairwise.complete.obs")
cormat %>% ggcorrplot(hc.order = TRUE, tl.cex = 5,
                      outline.col = "white",
                      insig = "blank") + ggtitle("Pearson correlation")


# correlation between mutations
first_half <- correlations_df_to_consider %>% 
  mutate("mut_cluster" = str_c(Gene_uniq2, cluster, sep = " ")) %>% 
  select("mutant" = Gene_uniq1, mut_cluster, correlation)
second_half <- correlations_df_to_consider %>% 
  mutate("mut_cluster" = str_c(Gene_uniq1, cluster, sep = " ")) %>% 
  select("mutant" = Gene_uniq2, mut_cluster, correlation)
mut_corr_df <- bind_rows(first_half, second_half)
corr_spread <- mut_corr_df %>% 
  spread(mutant, correlation)
corr_spread_matrix <- as.matrix(corr_spread[, -1])
rownames(corr_spread_matrix) <- corr_spread$mut_cluster
##### correlations between mutants
cormat <- cor(corr_spread_matrix[,], use = "pairwise.complete.obs")
cormat %>% ggcorrplot(hc.order = TRUE, tl.cex = 5,
                      outline.col = "white",
                      insig = "blank") + ggtitle("Pearson correlation")

representative_corr_df <- correlations_df %>% 
  filter(Gene_uniq1 %in% selected_mutants & Gene_uniq2 %in% selected_mutants) %>% 
  filter(cluster %in% representative_clusters$representative_clusters)
first_half <- representative_corr_df %>% 
  mutate("mut_cluster" = str_c(Gene_uniq2, cluster, sep = " - ")) %>% 
  select("mutant" = Gene_uniq1, mut_cluster, correlation)
second_half <- representative_corr_df %>% 
  mutate("mut_cluster" = str_c(Gene_uniq1, cluster, sep = " - ")) %>% 
  select("mutant" = Gene_uniq2, mut_cluster, correlation)
representative_mut_corr_df <- bind_rows(first_half, second_half)
representative_spread <- representative_mut_corr_df %>% 
  select(mutant, mut_cluster, correlation) %>% 
  spread(mut_cluster, correlation)
representative_spread_matrix <- as.matrix(representative_spread[-1])
rownames(representative_spread_matrix) <- representative_spread$mutant

pca_wrap <- function(data, pcn = 2) {
  #### first remove columns with more than half NA
  #data <- data[, -which(colMeans(is.na(data)) > 0.5)]  # remove columns that are more than half NA
  data <- data[, colSums(is.na(data)) < (nrow(data)/2)]
  data_preped <- pcaMethods::prep(data[, -1], center = FALSE, scale = "none")
  pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  scores <- as_tibble(scale(scores(pca))) %>% 
    add_column(., "mutant" = data$mutant)
  loadings <- as_tibble(scale(loadings(pca)))
  loadings <- add_column(loadings, "variables" = colnames(data[-1]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings, "cum_var" = R2cum(pca)))
}
principal_components <- str_c("PC", 1:10)
core_pca <- pca_wrap(representative_spread, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  separate(col = variables, into = c("mut", "cluster"), sep = " - ")
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "cluster"), 
               size = 5, alpha = 0.5) +
    geom_label(aes(label = mutant))
}
pdf("mutant_pairwise_corr_PCA.pdf", width = 15)
print(plots)
dev.off()
