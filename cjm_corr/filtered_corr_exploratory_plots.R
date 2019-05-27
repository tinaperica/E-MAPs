
library(tidyverse)

load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')


pdf('cjm_corr/hist_pval_corrs222.pdf')
x <- filtered_correlations %>% 
  select(query_uniq1, query_uniq2, pearson, two_sided_p_value) %>% 
  mutate('sig' = case_when((two_sided_p_value > 0.05) ~ FALSE, TRUE ~ TRUE))
ggplot(x, aes(x=pearson)) +
  geom_histogram(data = filter(x, sig == F), fill = 'red', alpha = 0.2, bins = 5000) +
  geom_histogram(data = filter(x, sig == T), fill = 'blue', alpha = 0.2, bins = 5000) +
  xlim(c(-0.25, 0.25)) + ggtitle('p-values of correlations of S-scores')
dev.off()

pdf('cjm_corr/srm1.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','srm1-ts')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `srm1-ts`)) + geom_point()
dev.off()

pdf('cjm_corr/srm1_rna1.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','rna1-1')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `rna1-1`)) + geom_point()
dev.off()

pdf('cjm_corr/srm1_acf4.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','acf4')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `acf4`)) + geom_point()
dev.off()



