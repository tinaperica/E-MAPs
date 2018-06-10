library(tidyverse)
library(mclust)
##### correlations between mutants
load("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/clustered_correlations/correlation_RData/mutants_only.RData")
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","NTER3XFLAG,WT","CTER3XFLAG,WT","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T")
cluster_annotation <- read_tsv("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/clustered_correlations/clusters/GO_slims_2018-03-05_pearson_complete_clusters.txt", col_names = T)
mut_correlations_df <- as_tibble(correlations_df) %>%
  filter(Gene_uniq1 %in% mutants & Gene_uniq2 %in% mutants) %>%
  mutate(correlation = as.numeric(correlation), random_correlation = as.numeric(random_correlation),
         random_correlation_2 = as.numeric(random_correlation_2), random_correlation_3 = as.numeric(random_correlation_3),
         random_correlation_4 = as.numeric(random_correlation_4), "pair" = str_c(Gene_uniq1, " ", Gene_uniq2))
pairs <- mut_correlations_df %>% pull(pair) %>% unique  
gathered_corr_data <- mut_correlations_df %>%  
  select(Gene_uniq1, Gene_uniq2, cluster, correlation, random_correlation, random_correlation_2, random_correlation_3, random_correlation_4) %>%
  gather(key = "correlation_type", value = "corr", -cluster, -Gene_uniq1, -Gene_uniq2)

corr_diff <- mut_correlations_df %>% 
  mutate("corr_diff1" = correlation - random_correlation, "corr_diff2" = correlation - random_correlation_2, 
         "corr_diff3" = correlation - random_correlation_3, "corr_diff4" = correlation - random_correlation_4) %>%
  group_by(cluster) %>%
  summarise(
    "correlation_whisker_diff" = boxplot(correlation)$stats[5, ] - boxplot(correlation)$stats[1, ], 
    "random_correlation_whisker_diff" = boxplot(c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)
    )$stats[5, ] - 
      boxplot( c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)
      )$stats[1, ],
    "correlation_upper_whisker" = boxplot(correlation)$stats[5, ],
    "random_upper_whisker" = boxplot(c(random_correlation, random_correlation_2, random_correlation_3, 
                                       random_correlation_4))$stats[5,],
    "diff_of_whisker_diffs" = correlation_whisker_diff - random_correlation_whisker_diff,
    "diff_of_upper_whiskers" = correlation_upper_whisker - random_upper_whisker,
    "correlation_sd" = sd(correlation),
    "random_sd" = sd(c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)),
    "sd_diff" = correlation_sd - random_sd
  ) %>%
  arrange(desc(sd_diff))
clusters.ordered.by.significance <- corr_diff %>% pull(cluster) %>% unique()
distribution_differences <- mut_correlations_df %>%
  group_by(cluster) %>%
  summarize("pvalue" = t.test(correlation, random_correlation)$p.value)
corr_diff <- inner_join(corr_diff, distribution_differences, by = "cluster")

corr_diff %>%
  select(cluster, sd_diff, diff_of_upper_whiskers, diff_of_whisker_diffs, pvalue) %>%
  ggplot(aes(x = sd_diff, y = diff_of_upper_whiskers, color = pvalue)) +
  geom_point() + geom_hline(yintercept = 0.4) + geom_vline(xintercept = 0.2) + 
  scale_colour_gradientn(colours = terrain.colors(10))

arranged_by_cluster_sig <- gathered_corr_data %>%
  mutate("category" = factor(cluster, levels = clusters.ordered.by.significance)) %>%
  arrange(category)
#plot <- arranged_by_cluster_sig %>%  ggplot(aes(x = corr, color = correlation_type)) + geom_bar() + facet_wrap(~category)
#ggsave("choose_clusters/many_random_vs_real_Gsp1_mut_correlations_by_library_cluster_by_sd.pdf", width = 48, height = 20)

clusters_to_consider <- corr_diff %>%
  filter(sd_diff > 0.2 & diff_of_upper_whiskers > 0.4 & pvalue > 0.25) %>%
  pull(cluster)
write_tsv(data.frame(clusters_to_consider), path = "~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/clusters_to_consider.txt")
arranged_by_cluster_sig_selected <- arranged_by_cluster_sig %>%
  filter(cluster %in% clusters_to_consider)
plot <- arranged_by_cluster_sig_selected %>% ggplot(aes(x = corr, color = correlation_type), alpha = 0.1) +
  geom_density() + facet_wrap(~category)
ggsave("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/many_random_vs_real_Gsp1_mut_correlations_SELECTED_clusters.pdf", width = 48, height = 20)

plot <- arranged_by_cluster_sig_selected %>% ggplot(aes(x = corr, color = correlation_type), alpha = 0.1) +
  geom_bar() + facet_wrap(~category)
ggsave("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/many_random_vs_real_Gsp1_mut_correlations_SELECTED_clusters_barplots.pdf", width = 48, height = 20)


selected_cluster_annotation <- cluster_annotation %>%
  filter(cluster %in% clusters_to_consider)
write_tsv(selected_cluster_annotation, path = "~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/selected_cluster_annotations.txt")


selected_cluster_annotation %>% 
  group_by(cluster) %>%
  summarize("count" = n()) %>%
  ggplot(aes(x = count)) + geom_histogram() + ggtitle("genes per cluster")
selected_cluster_annotation %>%
  group_by(ORF) %>%
  summarise("count" = n()) %>%
  ggplot(aes(x = count)) + geom_histogram() + ggtitle("clusters per gene")
genes_in_clusters <- cluster_annotation %>% pull(gene_name) %>% unique()
length(genes_in_clusters)
genes_in_selected_clusters <- selected_cluster_annotation %>% pull(gene_name) %>% unique()
length(genes_in_selected_clusters)
