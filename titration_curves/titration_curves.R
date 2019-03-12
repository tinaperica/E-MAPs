library(tidyverse)
library(cowplot)
library(viridis)
### complex annotation  - cluster the point mutations
load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
e.map <- as_tibble(e.map) %>% 
  select(-library_gene_name) %>% 
  mutate("mutant" = as.character(mutant), "library_ORF" = as.character(library_ORF)) %>% 
  mutate("mutant" = ifelse(mutant == "GSP1-NAT", "WT", mutant)) %>% 
  filter(! mutant %in% c("NTER3XFLAG_WT", "CTER3XFLAG_WT"))

#### Wodak complexes
##### keep only complexes that have more than one member that has a significant E-MAP score with at least one mutant
library_ORFs <- e.map %>% 
  filter(score < -2 | score > 2) %>% 
  pull(library_ORF) %>% unique()
complex_tibble <- read_tsv("titration_curves/all_wodak_complexes.txt")
complexes <- complex_tibble %>% 
  filter(ORF %in% library_ORFs) %>% 
  group_by(Complex) %>% 
  summarize("count" = n()) %>% 
  filter(count > 1) %>% pull(Complex)
complex_tibble <- complex_tibble %>% 
  filter(Complex %in% complexes)
###### after this filtering ended up with 115 complexes
emap_complex_data <- complex_tibble %>% 
  filter(Complex %in% complexes) %>% 
  inner_join(e.map, by = c("ORF" = "library_ORF")) %>% 
  select(-ORF) %>% 
  arrange(Name, mutant)




##### GAP/GEF ratio
GAP_GEF_ratio <- read_tsv("titration_curves/GAP_GEF_ratio.txt")
GAP_GEF_ratio_mutants <- GAP_GEF_ratio %>% pull(mutant) %>% unique()

for (i in seq_along(complexes)) {
  complex <- complexes[i]
  temp.emap <- emap_complex_data %>% 
    filter(Complex == complex) %>% 
    select(-Complex) %>% 
    spread(Name, score)
  mutants <- temp.emap %>% pull(mutant) %>% unique()
  temp.emap.matrix <- as.matrix(temp.emap[-1])
  rownames(temp.emap.matrix) <- temp.emap$mutant
  temp.emap.matrix[is.na(temp.emap.matrix)] <- 0 
  distance <- dist(temp.emap.matrix, method = "euclidean")
  mds <- cmdscale(distance, eig = T, k = 2)
  mds.coord <- as.tibble(mds$point) %>% 
    mutate("mutant" = mutants) %>% 
    rename("dim1" = V1, "dim2" = V2)
  ordered_genes.dim1 <- mds.coord %>% arrange(dim1) %>% pull(mutant)
  ordered_genes.dim2 <- mds.coord %>% arrange(dim2) %>% pull(mutant)
  mds_plot <- mds.coord %>%  ggplot(aes(dim1, dim2, label = mutant)) + 
    geom_text(check_overlap = TRUE)
  mean_wt_score <- emap_complex_data %>% 
    filter(Complex == complex & mutant == "WT") %>%
    summarize("mean" = mean(score)) %>% pull(mean)
  e.map_plot <- temp.emap %>% 
    gather(key = "complex member", value = "E-MAP score", -mutant) %>% 
    mutate("Gsp1 mutant" = factor(mutant, ordered_genes.dim1)) %>% 
    ggplot(aes(x = `Gsp1 mutant`, y = `E-MAP score`, color = `complex member`)) + 
    geom_point(size = 4, alpha = 0.75) +
    geom_hline(yintercept = mean_wt_score) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") +
    scale_color_viridis(discrete = T) + ggtitle(complex)
  
  titration_plot <- plot_grid(mds_plot, e.map_plot, nrow = 2, align = "v", rel_heights = c(1, 2))
  save_plot(titration_plot, 
         filename = str_c("titration_curves/output/titration_curves_all_mutants/",
           str_replace_all(complex, pattern = "/", replacement = ""), ".pdf"), base_width = 10, base_height = 7)
  
  e.map.GAPGEF <- temp.emap %>% filter(mutant %in% GAP_GEF_ratio_mutants)
  e.map_plot <- e.map.GAPGEF %>% 
    gather(key = "complex member", value = "E-MAP score", -mutant) %>% 
    mutate("Gsp1 mutant" = factor(mutant, ordered_genes.dim1[ordered_genes.dim1 %in% GAP_GEF_ratio_mutants])) %>% 
    ggplot(aes(x = `Gsp1 mutant`, y = `E-MAP score`, color = `complex member`)) + 
    geom_point(size = 4, alpha = 0.75) +
    geom_point(data = GAP_GEF_ratio, aes(x = mutant, y = log(GAP_GEF_ratio)), color = "red", alpha = 0.9) +
    geom_hline(yintercept = mean_wt_score) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") +
    scale_color_viridis(discrete = T) + ggtitle(complex)
  save_plot(e.map_plot, 
         filename = str_c("titration_curves/output/titration_curves_GAP_GEF/",
        str_replace_all(complex, pattern = "/", replacement = ""), ".pdf"), base_width = 10, base_height = 7)
  e.map_plot <- e.map.GAPGEF %>% 
    gather(key = "complex member", value = "E-MAP score", -mutant) %>% 
    mutate("Gsp1 mutant" = factor(mutant, GAP_GEF_ratio_mutants)) %>% 
    ggplot(aes(x = `Gsp1 mutant`, y = `E-MAP score`, color = `complex member`)) + 
    geom_point(size = 4, alpha = 0.75) +
    geom_point(data = GAP_GEF_ratio, aes(x = mutant, y = log(GAP_GEF_ratio)), color = "red", alpha = 0.9) +
    geom_hline(yintercept = mean_wt_score) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") +
    scale_color_viridis(discrete = T) + ggtitle(complex)
  save_plot(e.map_plot, 
         filename = str_c("titration_curves/output/titration_curves_GAP_GEF_order/",
                          str_replace_all(complex, pattern = "/", replacement = ""), ".pdf"), base_width = 10, base_height = 7)
}


prot_corr <- function(gene1, gene2) {
  corr <- cor(temp$score[temp$Name == gene1], temp$score[temp$Name == gene2], use = "pairwise.complete.obs")
}
mean_complex_correlation <- tibble("complex" = complexes, "mean_corr" = vector("double", length(complexes)) )
for (i in seq_along(complexes)) {
  temp <- emap_complex_data %>% 
    filter(Complex == complexes[i])
  complex_members <- temp %>% pull(Name) %>% unique()
  complex_member_pairs <- combn(complex_members, m = 2)
  complex_member_pairs <- list("gene1" = complex_member_pairs[1,], "gene2" = complex_member_pairs[2,])
  mean_complex_correlation$mean_corr[[i]] <- complex_member_pairs %>% pmap(prot_corr) %>% unlist() %>% mean
}
mean_complex_correlation %>% arrange(mean_corr) %>% ggplot(aes(x = mean_corr)) + geom_density()

complexes_to_keep <- mean_complex_correlation %>% filter(mean_corr > 0.5) %>% pull(complex)

selected_emap <- emap_complex_data %>% 
  filter(Complex %in% complexes_to_keep & mutant %in% GAP_GEF_ratio_mutants) %>% 
  group_by(Complex, mutant) %>% 
  mutate("mean_score" = mean(score, na.rm = T)) %>% 
  mutate("mean_score" = ifelse(mean_score < -5, -5, mean_score)) %>% 
  mutate("mean_score" = ifelse( (mean_score > -3 & mean_score < 2), 0, mean_score)) %>% 
  ungroup() %>% 
  select(-Name, -score) %>% unique()
complexes_to_keep <- selected_emap %>% group_by(Complex) %>% 
  summarize("mean" = mean(mean_score, na.rm = T)) %>% 
  filter(mean != 0) %>% pull(Complex)
selected_emap <- selected_emap %>% filter(Complex %in% complexes_to_keep)

bin_df <- tibble("score_bin_breaks" = c(-30, -7, -3, 2, 10),
                 "score_bin" = c(1, 1, 2, 3, 4),
                 "bin" = c(-2, -2, -1, 0, 2))
binned_emap <- selected_emap %>% 
  mutate("score_bin" = .bincode(mean_score, breaks = bin_df$score_bin_breaks, FALSE, FALSE)) %>% 
  inner_join(., bin_df, by = "score_bin") %>% 
  select(-score_bin_breaks, -score_bin, -mean_score) %>% unique()



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

spread_complexes <- binned_emap %>% 
  spread(Complex, bin)
ordered_complexes <- get_order_by_corr(spread_complexes)
spread_mutants <- binned_emap %>% 
  spread(mutant, bin)
ordered_mutants <- get_order_by_corr(spread_mutants)
binned_emap %>% 
    mutate("mutant" = factor(mutant, ordered_mutants)) %>%
    mutate("Complex" = factor(Complex, ordered_complexes)) %>% 
    arrange(mutant, Complex) %>% 
    ggplot(aes(x = mutant, y = Complex, fill = bin)) +
    geom_tile() + scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

spread_complexes <- selected_emap %>% 
  spread(Complex, mean_score)
ordered_complexes <- get_order_by_corr(spread_complexes)
spread_mutants <- selected_emap %>% 
  spread(mutant, mean_score)
ordered_mutants <- get_order_by_corr(spread_mutants)
selected_emap %>% 
  mutate("mutant" = factor(mutant, ordered_mutants)) %>%
  mutate("Complex" = factor(Complex, ordered_complexes)) %>% 
  arrange(mutant, Complex) %>% 
  ggplot(aes(x = mutant, y = Complex, fill = mean_score)) +
  geom_tile() + scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

