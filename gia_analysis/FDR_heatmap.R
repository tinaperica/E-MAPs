### gia FDR heatmap

library(tidyverse)
data <- read_tsv("gia_analysis/avg_merged_June2016_screen_for_Gia_gia_fdr.txt") %>% 
  select("mutant" = `Group A`, "gene_set" = `Group B`, FDR)

sig_gene_sets <- data %>% pull(gene_set) %>% unique()

gene_set <- read_tsv("gia_analysis/gia_right_annotations_groups.txt") %>% 
  filter(term %in% sig_gene_sets)
write_tsv(gene_set, "gia_analysis/gia_right_annotations_groups_preselected.txt")
data <- data %>% 
  group_by(gene_set) %>% 
  mutate("min_FDR" = min(FDR)) %>% 
  ungroup() %>%
  filter(min_FDR < 0.05) %>% 
  group_by(mutant) %>% 
  mutate("min_FDR" = min(FDR)) %>% 
  ungroup() %>% 
  filter(min_FDR < 0.05)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
get_reordered_cormat <- function(data) {
  cormat <- round(cor(data[, -1], use = "pairwise.complete.obs", method = "spearman"), 2)
  cormat <- reorder_cormat(cormat)
  return(cormat)
}
get_order_by_corr <- function(data) {
  cormat <- get_reordered_cormat(data)
  ordered <- rownames(cormat)
}

ordered_mut <- data %>% 
  spread(mutant, FDR) %>% 
  get_order_by_corr()
ordered_gene_set <- data %>% 
  spread(gene_set, FDR) %>% 
  get_order_by_corr()

data %>% 
  mutate("mutant" = factor(mutant, ordered_mut), "gene_set" = factor(gene_set, ordered_gene_set)) %>% 
  arrange(mutant, gene_set) %>% 
  ggplot(aes(x = mutant, y = gene_set, fill = FDR)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "white", midpoint = 0.075) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gsp1 mutant") + ylab("functional gene groups")
ggsave("gia_analysis/mutants_gene_sets_Sscore_enrichment_FDR_heatmap.pdf", height = 10, width = 10)

