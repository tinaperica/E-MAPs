library(tidyverse)
library(broom)
library(factoextra)
library(gplots)
library(pcaMethods)
library(GGally)
library(ggcorrplot)
library(fastcluster)

### correlation matrix for core subset
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}
gene_orf <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F) %>% 
  select("ORF" = X1, "gene_name" = X2) %>% 
  unique()
Gene_uniq_matcher <- read_tsv("basic_E-MAP_data/Gene_uniq_matcher.txt")
corr_of_corr_path <- "~/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Sep2018_analysis/corr_of_corr_mut_ms_partner_per_cluster/"
files <- dir(corr_of_corr_path)
(files <- files[grepl(files, pattern = ".RData")])
combined_data <- data.frame()
for (i in seq_along(files))  {    ## this loops through files (one per cluster)
  filepath <- file.path(corr_of_corr_path, files[i])
  load(filepath)
  corr_of_corr <- as_tibble(corr_of_corr) %>% 
    filter(corr_of_corr_type == "w_pearson_corr") %>% 
    mutate("fdr" = p.adjust(p.value, method = "fdr")) %>% 
    #filter(fdr < 0.05) %>% 
    inner_join(., Gene_uniq_matcher, by = c("geneB" = "Gene_uniq")) %>% 
    select(geneA, geneB, Gene.gene_name, cluster, corr, fdr)
  combined_data <- bind_rows(combined_data, corr_of_corr)
}
combined_data <- as_tibble(combined_data) %>% 
  rename("geneB_uniq" = geneB, "geneB" = Gene.gene_name) %>% 
  mutate("type" = ifelse(grepl(geneB_uniq, pattern = "_"), gsub(geneB_uniq, pattern = "^.+_", replacement = "", perl = T), "")) %>% 
  mutate("geneB_merge" = ifelse(type != "", str_c(geneB, type, sep = "_"), geneB)) %>% 
  select(-type)

save(combined_data, file = "reverse_corr_of_corr/20180925_combined_corr_of_corr.RData")
 
load("reverse_corr_of_corr/20180925_combined_corr_of_corr.RData")

spread_geneB_and_cluster <- combined_data %>% 
  mutate("geneB_cluster" = str_c(geneB_uniq, cluster, sep = " ")) %>% 
  select(geneB_cluster, geneA, corr) %>% 
  spread(geneB_cluster, corr)
spread_geneB_and_cluster <- data.frame(spread_geneB_and_cluster)
rownames(spread_geneB_and_cluster) <- spread_geneB_and_cluster$geneA
spread_geneB_and_cluster <- spread_geneB_and_cluster[, -1]
spread_geneB_and_cluster <- Filter(function (x) !all (is.na(x)), spread_geneB_and_cluster)  # this removes columns that are all NA
head(spread_geneB_and_cluster[1:10, 1:10])
hclust_mut_data <- hcut(spread_geneB_and_cluster, k = 8)
fviz_dend(hclust_mut_data, main = "mutants clustered by all corr of corr, k = 8", cex = 0.6)  ## from ‘factoextra’ package 
order_mutants <- hclust_mut_data$labels[hclust_mut_data$order]
mutant_groups <- cutree(hclust_mut_data, h = 250)
group_2_ordered_mutants <- order_mutants[order_mutants %in% names(mutant_groups[mutant_groups == 2])]
spread_geneA_and_cluster <- combined_data %>% 
  mutate("geneA_cluster" = str_c(geneA, cluster, sep = " ")) %>% 
  select(geneA_cluster, geneB_merge, corr) %>% 
  spread(geneA_cluster, corr)
spread_geneA_and_cluster <- data.frame(spread_geneA_and_cluster)
rownames(spread_geneA_and_cluster) <- spread_geneA_and_cluster$geneB_merge
spread_geneA_and_cluster <- spread_geneA_and_cluster[, -1]
spread_geneA_and_cluster <- Filter(function (x) !all (is.na(x)), spread_geneA_and_cluster)  # this removes columns that are all NA
head(spread_geneA_and_cluster[1:10, 1:10])
hclust_data <- hcut(spread_geneA_and_cluster, k = 2)
fviz_dend(hclust_data, main = "MS hits clustered by all corr of corr, k = 2", cex = 0.4)  ## from ‘factoextra’ package 
order_partners <- hclust_data$labels[hclust_data$order]
partner_groups <- cutree(hclust_data, h = 180)
group_2_ordered_partners <- order_partners[order_partners %in% names(partner_groups[partner_groups == 2])]
 

clusters <- combined_data %>% 
  pull(cluster) %>% unique()
partners_to_color <- combined_data %>% 
  filter(geneB %in% c("RNA1", "SRM1")) %>% 
  pull(geneB_merge) %>% unique()
coloring_list <- list()
colors <- c( "darkviolet", "darkcyan")
for ( i in seq_along(partners_to_color) ) {
  coloring_list[partners_to_color[i]] <- colors[i]
}
pdf("reverse_corr_of_corr/20180926_reverse_corr_of_corr_w_pearson_heatmaps_APMS_hits.pdf", height = 28, width = 14)
for (cn in seq_along(clusters)) {
  clust <- clusters[cn]
  temp <- combined_data %>% 
    filter(cluster == clust) %>% 
    select(geneB_merge, corr, geneA) %>% 
    mutate("geneA" = factor(geneA, order_mutants)) %>% 
    #mutate("geneB_merge" = factor(geneB_merge, order_partners)) %>% 
    arrange(geneA) %>% 
    spread(geneB_merge, corr)
  temp <- temp[, which(colMeans(is.na(temp)) < 0.5)] %>%   # remove columns that are more than half NA
    replace(is.na(.), 0)
  if (ncol(temp) > 6) {
    temp_df <- data.frame(temp)
    rownames(temp_df) <- temp_df$geneA
    temp.mat <- t(as.matrix(temp_df[, -1]))
    cols <- rep('black', nrow(temp.mat))
    for (i in seq_along(names(coloring_list))) {
      cols[row.names(temp.mat) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
    }
    heatmap.2(temp.mat, Rowv = F, Colv = F, col = cm.colors, scale = "none", trace = "none", cexCol = 0.8, 
              cexRow = 0.8, keysize = 0.5, margin = c(10, 10), colRow = cols, main = clust)
  }
}
dev.off()

pdf("reverse_corr_of_corr/20180926_reverse_corr_of_corr_w_pearson_heatmaps_GAP_GEF_groups.pdf", height = 10, width = 10)
for (cn in seq_along(clusters)) {
  clust <- clusters[cn]
  temp <- combined_data %>% 
    filter(cluster == clust & geneB_merge %in% group_2_ordered_partners) %>% 
    select(geneB_merge, corr, geneA) %>% 
    mutate("geneA" = factor(geneA, order_mutants)) %>% 
    #mutate("geneB_merge" = factor(geneB_merge, group_2_ordered_partners)) %>% 
    arrange(geneA) %>% 
    spread(geneB_merge, corr)
  temp <- temp[, which(colMeans(is.na(temp)) < 0.5)] %>%   # remove columns that are more than half NA
    replace(is.na(.), 0)
  if (ncol(temp) > 6) {
    temp_df <- data.frame(temp)
    rownames(temp_df) <- temp_df$geneA
    temp.mat <- t(as.matrix(temp_df[, -1]))
    cols <- rep('black', nrow(temp.mat))
    for (i in seq_along(names(coloring_list))) {
      cols[row.names(temp.mat) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
    }
    heatmap.2(temp.mat, Rowv = T, Colv = F, col = cm.colors, scale = "none", trace = "none", cexCol = 0.8, 
              cexRow = 0.8, keysize = 0.5, margin = c(10, 10), colRow = cols, main = clust)
  }
}
dev.off()

###
pca_wrap <- function(data, pcn = 2) {
  #### first remove columns with more than half NA
  data <- data[, which(colMeans(is.na(data)) < 0.2)]  # remove columns that are more than 20% NA
  data_preped <- pcaMethods::prep(data[, -1], center = FALSE, scale = "none")
  pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  scores <- as_tibble(scale(scores(pca))) %>% 
    add_column(., "mutant" = data$geneA)
  loadings <- as_tibble(scale(loadings(pca)))
  loadings <- add_column(loadings, "variables" = colnames(data[-1]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings, "cum_var" = R2cum(pca)))
}

core_regulators <- c("SRM1", "RNA1")
core_regulator_data <- combined_data %>% 
  filter(geneB %in% core_regulators) %>% 
  filter(geneA %in% group_2_ordered_mutants) %>% 
  mutate("cluster_gene_combo" = str_c(cluster, geneB_merge, sep = " ")) %>% 
  select(geneA, corr, cluster_gene_combo) %>% 
  spread(cluster_gene_combo, corr)

principal_components <- str_c("PC", 1:7)
core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  mutate("gene" = str_sub(variables, -11))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    geom_label(aes(label = mutant), label.size = 0.01) +
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "gene"), 
               size = 3, alpha = 0.2)
}
pdf("reverse_corr_of_corr/20180926_GAP_GEF_NIPALS_PCA.pdf")
print(plots)
dev.off()


core_regulators <- c("SRM1", "RNA1", "MOG1", "YRB1")
core_regulator_data <- combined_data %>% 
  filter(geneB %in% core_regulators) %>% 
  filter(geneA %in% group_2_ordered_mutants) %>% 
  mutate("cluster_gene_combo" = str_c(cluster, geneB, sep = " ")) %>% 
  select(geneA, corr, cluster_gene_combo) %>% 
  spread(cluster_gene_combo, corr)

principal_components <- str_c("PC", 1:7)
core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  mutate("gene" = str_sub(variables, -4))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    geom_label(aes(label = mutant), label.size = 0.01) +
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "gene"), 
               size = 3, alpha = 0.2)
}
pdf("reverse_corr_of_corr/20180926_GAP_GEF_MOG1_YRB1_NIPALS_PCA.pdf")
print(plots)
dev.off()


core_regulators <- c("SRM1", "MOG1", "YRB1")
core_regulator_data <- combined_data %>% 
  filter(geneB %in% core_regulators) %>% 
  filter(geneA %in% group_2_ordered_mutants) %>% 
  mutate("cluster_gene_combo" = str_c(cluster, geneB, sep = " ")) %>% 
  select(geneA, corr, cluster_gene_combo) %>% 
  spread(cluster_gene_combo, corr)

principal_components <- str_c("PC", 1:7)
core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  mutate("gene" = str_sub(variables, -4))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    geom_label(aes(label = mutant), label.size = 0.01) +
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "gene"), 
               size = 3, alpha = 0.2)
}
pdf("reverse_corr_of_corr/20180926_GEF_MOG1_YRB1_NIPALS_PCA.pdf")
print(plots)
dev.off()

