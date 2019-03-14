library(tidyverse)
library(RColorBrewer)
colors <- brewer.pal(n = 6, name = "Set1")
load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
e.map <- as_tibble(e.map) %>% 
  select(-library_gene_name) %>% 
  mutate("mutant" = as.character(mutant), "library_ORF" = as.character(library_ORF)) %>% 
  mutate("mutant" = ifelse(mutant == "GSP1-NAT", "WT", mutant)) 
mutants <- e.map %>% pull(mutant) %>% unique()
pearson_corr <- function(mut1, mut2) {
  e.map.mut1 <- e.map[e.map$mutant == mut1,]
  e.map.mut2 <- e.map[e.map$mutant == mut2,]
  temp <- inner_join(e.map.mut1, e.map.mut2, by = "library_ORF")
  corr <- cor(temp$score.x, temp$score.y,
              use = "pairwise.complete.obs")
}
rank_corr <- function(mut1, mut2) {
  e.map.mut1 <- e.map[e.map$mutant == mut1,]
  e.map.mut2 <- e.map[e.map$mutant == mut2,]
  temp <- inner_join(e.map.mut1, e.map.mut2, by = "library_ORF")
  corr <- cor(temp$score.x, temp$score.y,
              use = "pairwise.complete.obs", method = "kendall")
}
pairs <- t(combn(mutants, 2))
pairs <- list("mut1" = as.character(pairs[,1]), "mut2" = as.character(pairs[,2]))
mutant_pearson_correlation <- tibble(
  "mut1" = pairs$mut1, 
  "mut2" = pairs$mut2,
  "correlation" = vector("double", length(pairs$mut1))
)
mutant_pearson_correlation$correlation <- pairs %>% pmap(pearson_corr) %>% unlist()

mutant_pearson_correlation <- mutant_pearson_correlation %>% 
  mutate("method" = "pearson")


mutant_kendall_correlation <- tibble(
  "mut1" = pairs$mut1, 
  "mut2" = pairs$mut2,
  "correlation" = vector("double", length(pairs$mut1))
)
mutant_kendall_correlation$correlation <- pairs %>% pmap(rank_corr) %>% unlist()

mutant_kendall_correlation <- mutant_kendall_correlation %>% 
  mutate("method" = "kendall")

mutant_corr <- bind_rows(mutant_kendall_correlation, mutant_pearson_correlation) %>% 
  arrange(desc(correlation))
temp <- mutant_corr %>% 
  mutate("muta" = mut2, "mutb" = mut1) %>% 
  select("mut1" = muta, "mut2" = mutb, correlation, method)
mutant_corr <- bind_rows(mutant_corr, temp) %>% 
  arrange(mut1, mut2)
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

spread_mut <- mutant_corr %>% 
  filter(method == "pearson") %>% 
  select(-method) %>% 
  spread(mut2, correlation)
ordered_mut <- get_order_by_corr(spread_mut)



mutant_corr %>% 
  filter(method == "pearson") %>% 
  mutate("mut1" = factor(mut1, ordered_mut), "mut2" = factor(mut2, ordered_mut)) %>% 
  arrange(mut1, mut2) %>% 
  ggplot(aes(x = mut1, y = mut2, fill = correlation)) + 
  geom_tile() + scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gsp1 mutant") + ylab("Gsp1 mutant")
ggsave("plot_pairwise_scatterplots/corr_between_all_mutants_heatmap.pdf", height = 10, width = 10)


mutant_corr %>% 
  filter(method == "kendall") %>% 
  mutate("mut1" = factor(mut1, ordered_mut), "mut2" = factor(mut2, ordered_mut)) %>% 
  arrange(mut1, mut2) %>% 
  ggplot(aes(x = mut1, y = mut2, fill = correlation)) + 
  geom_tile() + scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gsp1 mutant") + ylab("Gsp1 mutant")
ggsave("plot_pairwise_scatterplots/rank_corr_between_all_mutants_heatmap.pdf", height = 10, width = 10)

mutant_pearson_correlation <- mutant_pearson_correlation %>% 
  arrange(desc(correlation))
selected_groups <- c("anaphase-promoting complex", "nuclear pore complex", "tRNA processing", "Swr1p complex")
gene_groups <- read_tsv("GSEA_like_analysis/gene_groups.txt") %>% 
  filter(term %in% selected_groups) %>% 
  select("library_ORF" = "ORF", term)
xy.limits <- c(min(e.map$score, na.rm = T), max(e.map$score, na.rm = T))
plots <- list()
#for (i in 1:10) {
for (i in seq_along(mutant_pearson_correlation$mut1)) {
  mutant1 <- mutant_pearson_correlation$mut1[i]
  mutant2 <- mutant_pearson_correlation$mut2[i]
  temp <- e.map %>% 
    filter(mutant == mutant1)
  to_plot <- e.map %>% 
    filter(mutant == mutant2) %>% 
    left_join(., gene_groups, by = "library_ORF") %>% 
    inner_join(., temp, by = "library_ORF") %>% 
    mutate("term" = ifelse(is.na(term), "other", term))
  to_plot_other <- to_plot %>% filter(term == "other")
  to_plot_selected <- to_plot %>% filter(term != "other")
  plots[[i]] <- to_plot_other %>% ggplot(aes(x = score.x, y = score.y)) + 
    geom_point(alpha = 0.3) +
    geom_point(data = to_plot_selected, aes(color = term)) + 
    xlab(mutant1) + ylab(mutant2) + xlim(xy.limits) + ylim(xy.limits) +
    ggtitle(str_c("Pearson correlation = ", signif(mutant_pearson_correlation$correlation[i], digits = 3))) +
    labs(color = "gene set") + theme_classic() + scale_color_manual( values = colors)
}
pdf("plot_pairwise_scatterplots/output/all_pairwise_scatterplots.pdf", height = 5, width = 7)
print(plots)
dev.off()
