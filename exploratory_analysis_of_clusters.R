library(tidyverse)

##### correlations between mutants
load("clustered_correlations/test/20180307_1_correlations.RData")
mut_correlations_df <- as_tibble(correlations_df)
gathered_corr_data <- mut_correlations_df %>%  
  select(cluster, correlation, random_correlation) %>%
  gather(key = "correlation_type", value = "corr", -cluster)

corr_diff <- mut_correlations_df %>% 
  mutate("corr_diff" = correlation - random_correlation) %>%
  group_by(cluster) %>%
  summarise("median_corr_diff" = median(corr_diff), "mean_corr_diff" = mean(corr_diff), 
            "corr_diff_spread" = max(corr_diff) - min(corr_diff),
            "diff_of_corr_medians" = median(correlation) - median(random_correlation), 
            "diff_of_corr_means" = mean(correlation) - mean(random_correlation),
            "diff_of_corr_spreads" = (max(correlation) - min(correlation)) - (max(random_correlation) - min(random_correlation)),
            "correlation_whisker_diff" = boxplot(correlation)$stats[5, ] - boxplot(correlation)$stats[1, ], 
            "random_correlation_whisker_diff" = boxplot(random_correlation)$stats[5, ] - boxplot(random_correlation)$stats[1, ],
            "diff_of_whisker_diffs" = correlation_whisker_diff - random_correlation_whisker_diff
            ) %>%
  arrange(desc(diff_of_whisker_diffs))
clusters.ordered.by.significance <- corr_diff %>% pull(cluster) %>% unique()

corr_diff %>%
  select(cluster, corr_diff_spread, correlation_whisker_diff, diff_of_corr_spreads, diff_of_whisker_diffs, random_correlation_whisker_diff) %>%
  gather(key = "measure", value = "value", -cluster) %>%
  ggplot(aes(x = value)) + geom_histogram() + facet_wrap(~measure)

plot <- gathered_corr_data %>%
  mutate("category" = factor(cluster, levels = clusters.ordered.by.significance)) %>%
  arrange(category) %>%
  ggplot(aes(x = corr, color = correlation_type)) + geom_bar() + facet_wrap(~category)
ggsave("random_vs_real_Gsp1_mut_correlations_by_library_cluster.pdf", width = 48, height = 20)
##### the plots are ordered by difference upper and lower whisker difference (diff_of_whisker_diffs)

good_clusters <- corr_diff %>% 
  filter(diff_of_whisker_diffs > 0.75) %>% pull(cluster) %>% unique()
  
gathered_corr_data %>% filter(cluster %in% good_clusters) %>%
  mutate("category" = factor(cluster, levels = good_clusters)) %>%
  ggplot(aes(x = category, y = corr, color = correlation_type)) + 
  geom_boxplot() + coord_flip()

ggsave("random_vs_real_Gsp1_mut_correlations_boxplots.pdf", height = 45, width = 10)

