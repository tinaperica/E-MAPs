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
#load("20180507_corr_of_corr_all.RData")  # combined_data
corr_of_corr_path <- "clustered_correlations/corr_of_corr/"
files <- dir(corr_of_corr_path)
(files <- files[grepl(files, pattern = "soft_cos_sim_20180914_all_corr_of_corr.RData")])
combined_filtered_data <- data.frame()
for (i in seq_along(files))  {    ## this loops through files (one per cluster)
  filepath <- file.path(corr_of_corr_path, files[i])
  load(filepath)
  combined_filtered_data <- bind_rows(combined_filtered_data, soft_cos_sim)
}
rm(soft_cos_sim)

#load("20180507_corr_of_corr_filter_FDR.RData") ## combined_filtered_data
combined_filtered_data <- as_tibble(combined_filtered_data) %>% 
  separate(col = geneB, into = c("ORF"), sep = "_", remove = F) %>% 
  inner_join(., gene_orf, by = "ORF")
combined_filtered_data <- combined_filtered_data %>% 
  select(geneA, geneB, cluster, ORF, "corr" = corr.x, gene_name) %>% 
  mutate("corr" = as.numeric(corr))
  

### filtered data is filtered by FDR (filter(FDR < 0.05 & random_avrg_FDR > 0.05))
#### only keep cluster + geneB combinations where at least one of the mutants has a high correlation
### i.e. absolute value larger than 0.3
### but don't keep if the abs(corr) in WT negative control (GSP1-NAT) is larger than 0.3
# gene_cluster_pairs_with_strong_WT_corr <- combined_filtered_data %>% 
#   filter(geneA == "GSP1-NAT" & abs(corr) > 0.3) %>% 
#   mutate("uniq_gene_name" = str_c(geneB, gene_name, sep = " ")) %>% 
#   select(geneA, uniq_gene_name, cluster)
# strong_corr_subset <- combined_filtered_data %>% 
#   group_by(geneB, cluster) %>%
#   mutate("per_mutant_max_corr" = max(corr)) %>% 
#   filter(per_mutant_max_corr > 0.3) %>% 
#   ungroup() %>% 
#   mutate("uniq_gene_name" = str_c(geneB, gene_name, sep = " "))
# 
# strong_corr_subset <- anti_join(strong_corr_subset, 
#                     gene_cluster_pairs_with_strong_WT_corr, by = c("uniq_gene_name", "cluster"))
# 
# ###### use the strong data in all the downstream analysis
# combined_filtered_data <- strong_corr_subset

ms_hits <- read_tsv("SAINT_MS_hits.txt", col_names = F) %>% 
  rename("ORF" = X1) %>% 
  unique()
add <- c("YGR218W", "YER110C", "YAR002W", "YGL241W", "YDR192C", "YDR335W", "YJR074W", 
      "YGL097W", "YMR235C", "YKL205W", "YHR200W", "YMR308C", "YDR002W")
ms_hits <- bind_rows(ms_hits, tibble("ORF" = add)) %>% 
  unique()
# # sanity checks
# combined_data %>% pull(gene_name) %>% unique()
# combined_data %>% filter(gene_name == "MOG1" & str_detect(string = cluster, pattern = "nucle")) %>% 
#   ggplot(aes(x = geneA, y = corr, color = cluster)) + 
#   geom_point() + coord_flip() #+ theme(legend.position="none")
# # mean and variance for each cluster
# mean_and_var_per_cluster <- combined_data %>% 
#   group_by(cluster) %>% 
#   summarize("mean_per_cluster" = mean(corr, na.rm = T), "variance_per_cluster" = var(corr, na.rm = T)) %>% 
#   gather(measure, value, - cluster) 
# mean_and_var_per_cluster %>% ggplot(aes(x = value)) + geom_density() + facet_wrap(~measure)
# clusters_with_low_values <- mean_and_var_per_cluster %>% 
#   filter(abs(value) < 0.0005) %>% 
#   arrange(abs(value))
# 
# combined_data %>% 
#   filter(cluster == "RNA catabolic process") %>% 
#   ggplot(aes(x = corr)) + geom_histogram()
# combined_data %>% 
#   filter(cluster == "cell cortex_GO_1") %>% 
#   ggplot(aes(x = corr)) + geom_histogram()
# combined_data %>% 
#   filter(cluster == "RNA catabolic process" & corr > 0.5 & geneA == "GSP1-NAT")
# 
# ### check both clusters and genes that strongly correlate with the negative control GSP1-NAT
# WT_data <- combined_data %>% 
#   filter(geneA == "GSP1-NAT" & abs(corr) > 0.5)
# WT_data_filtered <- combined_filtered_data %>% 
#   filter(geneA == "GSP1-NAT" & abs(corr) > 0.5)
# #####


core_partners <- c("YRB2", "YRB1", "KAP95", "NUP1", "NUP60", "KAP104", "SRM1", "LOS1",
                   "NTF2", "CSE1", "CRM1", "RNA1", "MSN5", "PSE1", "SRP1", "MOG1", "MTR10", "NMD5",
                   "NUP159", "NUP42", "KAP120", "KAP123", "KAP114", "RPN10")
core_geneBs <- combined_filtered_data %>% 
  filter(gene_name %in% core_partners) %>% 
  pull(geneB) %>% unique()
core_subset <- combined_filtered_data %>% 
  filter(geneB %in% core_geneBs) %>% 
  mutate("geneB_cluster" = str_c(geneB, cluster, sep = " ")) %>% 
  select(geneB_cluster, geneA, corr) %>% 
  spread(geneB_cluster, corr)
core_subset <- Filter(function (x) !all (is.na(x)), core_subset)  # this removes columns that are all NA
core_subset <- data.frame(core_subset)
rownames(core_subset) <- core_subset$geneA
core_subset <- core_subset[, -1]
core_hclust_data <- hcut(dist(core_subset), k = 8)
fviz_dend(core_hclust_data, main = "mutants clustered by corr of corr with core partners k = 8", cex = 0.6)  ## from ‘factoextra’ package 


spread_geneB_and_cluster <- combined_filtered_data %>% 
  mutate("geneB_cluster" = str_c(geneB, cluster, sep = " ")) %>% 
  select(geneB_cluster, geneA, corr) %>% 
  spread(geneB_cluster, corr)
spread_geneB_and_cluster <- data.frame(spread_geneB_and_cluster)
rownames(spread_geneB_and_cluster) <- spread_geneB_and_cluster$geneA
spread_geneB_and_cluster <- spread_geneB_and_cluster[, -1]
spread_geneB_and_cluster <- Filter(function (x) !all (is.na(x)), spread_geneB_and_cluster)  # this removes columns that are all NA
head(spread_geneB_and_cluster[1:10, 1:10])
hclust_data <- hcut(spread_geneB_and_cluster, k = 8)
fviz_dend(hclust_data, main = "mutants clustered by all corr of corr, k = 8", cex = 0.6)  ## from ‘factoextra’ package 

order_of_mutants_based_on_all_corr_of_corr <- hclust_data$labels[hclust_data$order]
groups <- cutree(hclust_data, h = 500)
mut_in_group_with_WT <- names(groups[groups == 1])
spread_geneB_and_cluster_strong_mutants <- spread_geneB_and_cluster[
  !(row.names(spread_geneB_and_cluster) %in% mut_in_group_with_WT), ]
mutants_diff_from_WT <- names(groups[groups != 1])
hclust_data_strong <- hcut(spread_geneB_and_cluster_strong_mutants, k = 5)
fviz_dend(hclust_data_strong, main = "Strong mutants clustered by all the correlations of correlations", cex = 0.6)  ## from ‘factoextra’ package 


spread_geneA_and_cluster <- combined_filtered_data %>% 
  filter(ORF %in% ms_hits$ORF) %>% 
  mutate("gene_name_uniq" = str_c(gene_name, geneB, sep = " ")) %>% 
  mutate("geneA_cluster" = str_c(geneA, cluster, sep = " ")) %>% 
  select(geneA_cluster, gene_name_uniq, corr) %>% 
  spread(geneA_cluster, corr)
spread_geneA_and_cluster <- data.frame(spread_geneA_and_cluster)
rownames(spread_geneA_and_cluster) <- spread_geneA_and_cluster$gene_name_uniq
spread_geneA_cluster <- spread_geneA_and_cluster[, -1]
spread_geneA_and_cluster <- spread_geneA_and_cluster[, which(colMeans(is.na(spread_geneA_and_cluster)) < 0.5)]
spread_geneA_and_cluster <- spread_geneA_and_cluster[which(rowMeans(is.na(spread_geneA_and_cluster)) < 0.5), ]

spread_geneA_and_cluster <- Filter(function (x) !all (is.na(x)), spread_geneA_and_cluster)
head(spread_geneA_and_cluster[1:10, 1:10])
hclust_data <- hcut(spread_geneA_and_cluster, k = 6)
order_of_ms_hits_genes <- hclust_data$labels[hclust_data$order]
fviz_dend(hclust_data, main = "MS hits clustered by corr of corr with Gsp1 mutants, k = 10", cex = 0.5)  ## from ‘factoextra’ package 
ggsave("MS_hits_dendrogram_by_EMAP_corr_of_corr.pdf", width = 20)

groups <- cutree(hclust_data, h = 75)
GEF_group <- names(groups[groups == 7])
GAP_group <- names(groups[groups == 6])
GAP_and_GEF_groups <- names(groups[groups == 6 | groups == 7])
spread_geneB_and_cluster_strong_mutants <- spread_geneB_and_cluster[
  !(row.names(spread_geneB_and_cluster) %in% mut_in_group_with_WT), ]
mutants_diff_from_WT <- names(groups[groups != 1])
hclust_data_strong <- hcut(spread_geneB_and_cluster_strong_mutants, k = 5)
fviz_dend(hclust_data_strong, main = "Strong mutants clustered by all the correlations of correlations", cex = 0.6)  ## from ‘factoextra’ package 


cormat <- round(cor(t(core_subset), use = "pairwise.complete.obs"),2)
cormat <- reorder_cormat(cormat)
cormat %>% ggcorrplot(outline.col = "white",insig = "blank", tl.cex = 7) + 
  ggtitle("Pearson correlation, by core partners", "strong correlations subset") 


cormat <- round(cor(t(spread_geneB_and_cluster), use = "pairwise.complete.obs"),2)
cormat <- reorder_cormat(cormat)
cormat %>% ggcorrplot(outline.col = "white",insig = "blank", tl.cex = 7) + 
  ggtitle("Pearson correlation, all", "strong correlations subset") 


### make per cluster heatmaps from combined_filtered_data using heatmap.2
# make two heatmaps - one for core partners only and one for everything pulled down with APMS (ms_hits)


ms_hits_corr_data <- inner_join(combined_filtered_data, ms_hits, by = "ORF") %>% 
  select(cluster, geneA, geneB, gene_name, corr) %>% 
  mutate("uniq_gene_name" = str_c(gene_name, geneB, sep = " "))
clusters <- ms_hits_corr_data %>% 
  pull(cluster) %>% unique()

partners_to_color <- ms_hits_corr_data %>% 
  filter(geneB %in% core_geneBs) %>% 
  mutate("uniq_gene_name" = gsub(uniq_gene_name, pattern = " ", replacement = ".")) %>% 
  pull(uniq_gene_name) %>% unique()
coloring_list <- list()
#colors <- c("dodgerblue3", "brown3", "darkcyan", "darkviolet", "darkorange3", "darkgoldenrod2")
#colors <- c("coral1", "dodgerblue3", "deeppink3", "dodgerblue3", "deeppink3", 
 #           "plum4",  "brown3", "darkcyan", "darkviolet", "red3", "palevioletred4", "darkorange3", "darkgoldenrod2")
colors <- c("coral1", "coral1", "darkcyan", "darkviolet", 
           "coral1",  "coral1", "coral1", "coral1", "coral1", "coral1", "coral1", "red3", "coral1")

for ( i in seq_along(partners_to_color) ) {
  coloring_list[partners_to_color[i]] <- colors[i]
}
pdf("20180915_corr_of_corr_soft_cos_sim_heatmaps_APMS_hits.pdf", height = 25, width = 14)
for (cn in seq_along(clusters)) {
  clust <- clusters[cn]
  temp <- ms_hits_corr_data %>% 
    filter(cluster == clust) %>% 
    select(uniq_gene_name, corr, geneA) %>% 
    mutate("geneA" = factor(geneA, order_of_mutants_based_on_all_corr_of_corr)) %>% 
    #mutate("uniq_gene_name" = factor(uniq_gene_name, order_of_ms_hits_genes)) %>% 
    filter(! is.na(uniq_gene_name)) %>% 
    arrange(geneA, uniq_gene_name) %>% 
    spread(uniq_gene_name, corr)
  temp <- temp[, which(colMeans(is.na(temp)) < 0.8)] %>%   # remove columns that are more than half NA
    replace(is.na(.), 0)
  if (ncol(temp) > 6) {
    temp_df <- data.frame(temp)
    rownames(temp_df) <- temp_df$geneA
    temp.mat <- t(as.matrix(temp_df[, -1]))
    cols <- rep('black', nrow(temp.mat))
    for (i in seq_along(names(coloring_list))) {
      cols[row.names(temp.mat) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
    }
    heatmap.2(temp.mat, Colv = F, col = cm.colors, scale = "none", trace = "none", cexCol = 0.8, 
            cexRow = 0.8, keysize = 0.5, margin = c(10, 10), colRow = cols, main = clust)
  }
}
dev.off()

#### heatmaps with just core partners
pdf("20180915_corr_of_corr_soft_cos_sim_heatmaps_partners.pdf", height = 10, width = 15)
for (cn in seq_along(clusters)) {
  clust <- clusters[cn]
  temp <- ms_hits_corr_data %>%
    filter(cluster == clust & geneB %in% core_geneBs) %>% 
    select(uniq_gene_name, corr, geneA) %>% 
    spread(uniq_gene_name, corr) %>% 
    replace(is.na(.), 0)
  temp_df <- data.frame(temp)
  rownames(temp_df) <- temp_df$geneA
  temp.mat <- t(as.matrix(temp_df[, -1]))
  cols <- rep('black', nrow(temp.mat))
  for (i in seq_along(names(coloring_list))) {
    cols[row.names(temp.mat) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
  }
  heatmap.2(temp.mat, col = cm.colors, scale = "none", trace = "none", cexCol = 0.8, 
            cexRow = 0.8, key = F, margin = c(7, 10), colRow = cols, main = clust)
}
dev.off()

### gap and gef groups
#### heatmaps with just core partners
pdf("20180915_corr_of_corr_soft_cos_sim_heatmaps_GAP_and_GEF_group.pdf", height = 10, width = 15)
for (cn in seq_along(clusters)) {
  clust <- clusters[cn]
  temp <- ms_hits_corr_data %>%
    filter(cluster == clust & uniq_gene_name %in% GAP_and_GEF_groups) %>% 
    select(uniq_gene_name, corr, geneA) %>% 
    mutate("geneA" = factor(geneA, order_of_mutants_based_on_all_corr_of_corr)) %>% 
    spread(uniq_gene_name, corr) %>% 
    replace(is.na(.), 0)
  temp_df <- data.frame(temp)
  rownames(temp_df) <- temp_df$geneA
  temp.mat <- t(as.matrix(temp_df[, -1]))
  cols <- rep('black', nrow(temp.mat))
  for (i in seq_along(names(coloring_list))) {
    cols[row.names(temp.mat) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
  }
  heatmap.2(temp.mat, Colv = F, col = cm.colors, scale = "none", trace = "none", cexCol = 0.8, 
            cexRow = 0.8, key = F, margin = c(7, 10), colRow = cols, main = clust)
}
dev.off()


### core regulators
core_regulators <- c("MOG1 YJR074W", "YRB1 YDR002W_TSQ582", "RNA1 YMR235C_TSQ172", "SRM1 YGL097W_TSQ958")
ms_hits_corr_data %>% 
  filter(uniq_gene_name %in% core_regulators) %>% 
  mutate() %>% 
  mutate("geneB_cluster_combo" = str_c(gene_name, cluster, sep = " ")) %>% 
  select(geneA, geneB_cluster_combo, corr)
coloring_list <- list()
colors <- c("dodgerblue3", "darkgoldenrod2", "darkviolet", "darkcyan")

for ( i in seq_along(partners_to_color) ) {
  coloring_list[partners_to_color[i]] <- colors[i]
}
# for each mutant find the highest cluster/geneB combinations and explore those


## group clusters by correlation between clusters?

###
pca_wrap <- function(data, pcn = 2) {
  #### first remove columns with more than half NA
  data <- data[, which(colMeans(is.na(data)) < 0.5)]  # remove columns that are more than half NA
  data_preped <- pcaMethods::prep(data[, -1], center = FALSE, scale = "none")
  pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  scores <- as_tibble(scale(scores(pca))) %>% 
    add_column(., "mutant" = data$geneA)
  loadings <- as_tibble(scale(loadings(pca)))
  loadings <- add_column(loadings, "variables" = colnames(data[-1]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings, "cum_var" = R2cum(pca)))
}

# loadings <- loadings %>% 
#   mutate("gene" = str_sub(cluster_gene_pair, -4))
# scores_plots <- list()
# for (i in 2:pcn) {
#   PC_n <- as.character(i)
#   plots[[i]] <- scores %>% ggplot(., aes_string("PC1", PC_n), label = geneA) + geom_point() +
#     geom_label_repel(aes(label = mutant), box.padding   = 0.35, point.padding = 0.5,
#                      segment.color = 'grey50') + theme_classic()
#   loadings %>% ggplot(., aes_string(PC1, PC2, color = gene)) + geom_point() 
#   loadings %>% ggplot(., aes(PC1, PC3, color = gene)) + geom_point() 
#   loadings %>% ggplot(., aes(PC2, PC3, color = gene)) + geom_point() 
# }

# core_regulators <- c("YRB1", "KAP95", "SRM1", 
#                      "NTF2", "CSE1", "CRM1", "RNA1", "SRP1", "MOG1")
core_regulators <- c("SRM1", "RNA1", "MOG1", "YRB1")
core_regulator_geneBs <- combined_filtered_data %>% 
  filter(gene_name %in% core_regulators) %>% 
  pull(geneB) %>% unique()
core_regulator_data <- combined_filtered_data %>% 
  filter(geneB %in% core_regulator_geneBs) %>% 
  filter(! geneA %in% mut_in_group_with_WT) %>% 
  mutate("cluster_gene_combo" = str_c(cluster, gene_name, sep = " ")) %>% 
  select(geneA, corr, cluster_gene_combo) %>% 
  spread(cluster_gene_combo, corr)

sig_Gsp1_GI_15_7_mut_data <- ms_hits_corr_data %>% 
  filter(cluster == "sig_Gsp1_GI_15_7_mut") %>% 
  select(geneA, corr, uniq_gene_name) %>% 
  spread(uniq_gene_name, corr)

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
pdf("20180918_GAP_GEF_NIPALS_PCA_strong.pdf")
print(plots)
dev.off()


#### correlate the GO clusters
cluster_genes_spread <- combined_filtered_data %>% 
  mutate("gene_combo" = str_c(geneA, geneB, sep = " ")) %>% 
  select(gene_combo, cluster, corr) %>% 
  spread(gene_combo, corr)
cluster_genes_spread <- data.frame(cluster_genes_spread)
rownames(cluster_genes_spread) <- cluster_genes_spread$cluster
cluster_genes_spread <- cluster_genes_spread[,-1]
cormat <- round(cor(cluster_genes_spread, use = "pairwise.complete.obs"),2)
cormat <- reorder_cormat(cormat)
cormat %>% ggcorrplot(outline.col = "white",insig = "blank", tl.cex = 4) + 
  ggtitle("Go clusters", "Pearson correlation") 


hclust_data <- hcut(cluster_genes_spread, k = 8)
fviz_dend(hclust_data, main = "GO clusters hclustered by all corr of corr, k = 8", cex = 0.6)  ## from ‘factoextra’ package 

GO_cluster_groups <- data.frame("group" = cutree(hclust_data, h = 150))
GO_cluster_groups$cluster <- rownames(GO_cluster_groups)
GO_cluster_groups <- as_tibble(GO_cluster_groups) %>% 
  arrange(group)
# group_names <- tibble("group" = 1:13, "group name" = 
#               c("cell cycle", "cell_cortex_GO_1", 
#                 "DNA, transcription, chromatin, cell cycle", 
#                 "DNA repair, protein targeting", "chromatin, transcriptiom", 
#                 "chromatin binding", "chromatin modification", "DNA, lipids, mitochondria",
#                 "ribosomes, translation", "ER and Golgi", "ER", "histone modifications", 
#                 "nuclear transport"))
group_names <- tibble("group" = 1:27, "group name" = 
      c("cell cycle", "cell cortex", "cell cycle, chromosome", 
        "cell cycle, telomere", "DNA repair", "chromatin, transcription", 
        "chromatin binding", "peptidyl-amino acid modification", "DNA binding",
        "ribosomes, translation", "endomembrane system",  "protein targeting", 
        "ER", "ER", "histone modification", "lipids", "mitochondria", 
        "RNA processing", "nuclear transport", 
        "protein modification by small protein conjugation", "protein targeting", 
        "regulation of cell cycle and proteolysis", 
        "transcription and mRNA processing GO 6", "lib17", "lib19", "lib22", "lib24"))

GO_cluster_groups <- GO_cluster_groups %>% inner_join(., group_names, by = "group")

combined_filtered_data <- combined_filtered_data %>% inner_join(., GO_cluster_groups, by = "cluster")


all_spread <- combined_filtered_data %>% 
  mutate("gene_cluster_combo" = str_c(geneB, cluster, sep = " ")) %>% 
  mutate("combo" = str_c(gene_cluster_combo, `group name`, sep = " - ")) %>% 
  select(geneA, combo, corr) %>% 
  spread(combo, corr)
  
principal_components <- str_c("PC", 1:10)
core_pca <- pca_wrap(all_spread, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  separate(col = variables, into = c("combo", "cluster_group"), sep = " - ")
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "cluster_group"), 
               size = 5, alpha = 0.5) +
    geom_label(aes(label = mutant))
}
pdf("all_NIPALS_PCA_GO_cluster_groups.pdf", width = 12)
print(plots)
dev.off()



strong_mutants <- c("CTER3XFLAG WT","D79A","D79S",
                    "G80A","H141E","H141I",
                    "H141R","H141V","K101R",
                    "K143W",
                    "NTER3XFLAG WT","Q147E","R108A",
                    "R108G","R108I","R108L","R108Q","R108S","R108Y",
                    "R112A","R112S","R78K",
                    "T34A","T34E","T34G","T34L",
                    "T34Q","Y148I","Y157A")
select_groups <- c("nuclear transport", "cell cortex", "ER", "ribosomes, translation", "cell cycle",
                   "cell cycle, chromosome", 
                   "cell cycle, telomere", "DNA repair", "chromatin, transcription", 
                   "chromatin binding", "peptidyl-amino acid modification")
core_regulators <- c( "MOG1")
core_partners_spread <- combined_filtered_data %>% 
  filter(gene_name %in% core_regulators) %>% 
  mutate("gene_cluster_combo" = str_c(geneB, cluster, sep = " ")) %>% 
  mutate("combo" = str_c(`group name`, gene_name, gene_cluster_combo, sep = " - ")) %>% 
  select(geneA, combo, corr) %>% 
  spread(combo, corr)

principal_components <- str_c("PC", 1:10)
core_pca <- pca_wrap(core_partners_spread, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  separate(col = variables, into = c("cluster_group", "partner", "combo"), sep = " - ")
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "cluster_group", shape = "partner"), 
               size = 5, alpha = 0.5) +
    geom_label(aes(label = mutant))
}
pdf("MOG1_NIPALS_PCA_GO_cluster_groups.pdf", width = 12)
print(plots)
dev.off()






select_data <- combined_filtered_data %>% 
  filter(`group name` %in% select_groups) %>% 
  mutate("gene_cluster_combo" = str_c(geneB, cluster, sep = " ")) %>% 
  mutate("combo" = str_c(gene_cluster_combo, `group name`, sep = " - ")) %>% 
  select(geneA, combo, corr) %>% 
  spread(combo, corr)

core_pca <- pca_wrap(select_data, pcn = length(principal_components))
core_pca$loadings <- core_pca$loadings %>% 
  separate(col = variables, into = c("combo", "cluster_group"), sep = " - ")
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = core_pca$scores, aes_string("PC1", pc)) + 
    theme_classic() +
    geom_point(data = core_pca$loadings, aes_string("PC1", pc, color = "cluster_group"), 
               size = 5, alpha = 0.5) +
    geom_label(aes(label = mutant))
}
pdf("all_NIPALS_PCA_select_GO_cluster_groups.pdf", width = 12)
print(plots)
dev.off()

GO_cluster_group <- GO_cluster_groups %>% pull(`group name`) %>% unique()
for (i in seq_along(GO_cluster_group)) {
  group <- GO_cluster
}

ms_mutants <- c("A180T", "CTER3XFLAG WT", "D79A", "D79S", "F58A", "G80A", "H141E", "H141I", "H141R", "H141V", "K101R", "K143W",
            "NTER3XFLAG WT", "Q147E", "R108A", "R108G", "R108I", "R108L", "R108Q",    
            "R108Y", "R112A", "R112S", "R78K", "T34A", "T34E", "T34G", "T34L", "T34Q", "Y148I", "Y157A")

strong_mutants <- c("CTER3XFLAG WT","D79A","D79S",
                    "G80A","H141E","H141I",
                    "H141R","H141V","K101R",
                    "K143W",
                    "NTER3XFLAG WT","Q147E","R108A",
                    "R108G","R108I","R108L","R108Q","R108S","R108Y",
                    "R112A","R112S","R78K",
                    "T34A","T34E","T34G","T34L",
                    "T34Q","Y148I","Y157A")