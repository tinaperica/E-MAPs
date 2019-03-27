### gia FDR heatmap

library(tidyverse)
data <- read_tsv("gia_analysis/emap_screen_filtered_gia_fdr.txt") %>% 
  select("mutant" = `Group A`, "gene_set" = `Group B`, FDR) %>% 
  arrange(gene_set)

# gene_set_enriched_in_wt <- data %>% 
#   filter(mutant == "WT" & gene_set != "nuclear transport" & 
#   gene_set != "nuclear pore complex" & gene_set != "nucleobase-containing compound transport") %>%
#   pull(gene_set) %>% unique()
### [1] "nuclear pore complex"                     "transmembrane transport"                 
# [3] "nucleobase-containing compound transport" "mitochondrion organization"              
# [5] "peroxisome"                               "mitochondrion"                           
# [7] "nuclear transport"  

sig_gene_sets <- data %>% 
  filter(! gene_set %in% c('chromatin', 'chromosome', 'endomembrane system', 'histone binding', 
                      'invasive growth in response to glucose limitation', 'isomerase activity', 'kinase activity', 'mitochondria',
                      'protein targeting', 'response to heat')) %>% 
  pull(gene_set) %>% unique()

# sig_mut <- data %>% group_by(mutant) %>% 
#   mutate("nrich" = n()) %>% 
#   ungroup() %>% 
#   filter(nrich > 8) %>% pull(mutant) %>% unique()
gene_set <- read_tsv("gia_analysis/gia_right_annotations_gene_groups.txt") %>% 
  filter(term %in% sig_gene_sets)
write_tsv(gene_set, "gia_analysis/gia_right_annotations_groups_preselected.txt")
#data <- data %>% filter(mutant %in% sig_mut)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
get_reordered_cormat <- function(data) {
  cormat <- round(cor(data[, -1], use = "pairwise.complete.obs", method = "pearson"), 2)
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
  arrange(gene_set, mutant) %>% 
  ggplot(aes(x = gene_set, y = mutant, fill = FDR)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "white", midpoint = 0.075) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gsp1 mutant") + ylab("functional gene groups")

ggsave("gia_analysis/mutants_gene_sets_Sscore_enrichment_FDR_heatmap.pdf", height = 13, width = 12)

