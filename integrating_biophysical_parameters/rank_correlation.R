library(tidyverse)
library(viridis)
library(pracma)
library(RColorBrewer)
colors <- brewer.pal(n = 6, name = "Set1")
orf_gene_index <- read_tsv(file = "orf_gene_GO_sgd_annotation.txt", col_names = F) %>% 
  select("ORF" = X1, "gene_name" = X2) %>% unique()

bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')
##### bad strains are library strains to be removed (according to Amelie's file, and Hannes's instruction, those should be removed)

### complex annotation  - cluster the point mutations
load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
e.map <- as_tibble(e.map) %>% 
  select(-library_gene_name) %>%
  filter(! (library_ORF %in% bad_strains)) %>% 
  mutate("mutant" = as.character(mutant), "library_ORF" = as.character(library_ORF)) %>% 
  mutate("mutant" = ifelse(mutant == "GSP1-NAT", "WT", mutant)) %>% 
  filter(! mutant %in% c("NTER3XFLAG_WT", "CTER3XFLAG_WT"))

### only look at library genes that have a siginificant S-score score with at least one of the mutants
library_genes <- e.map %>% 
  filter(score < -2 | score > 2) %>% 
  pull(library_ORF) %>% unique()
e.map <- e.map %>% 
  filter(library_ORF %in% library_genes)

##### GAP/GEF ratio
GAP_GEF_ratio <- read_tsv("titration_curves/GAP_GEF_ratio_all.txt")
GAP_GEF_ratio_mutants <- GAP_GEF_ratio %>% pull(mutant) %>% unique()
GAP_GEF_ratio <- GAP_GEF_ratio %>% 
  gather(param, value, -mutant) %>% 
  group_by(param) %>% 
  mutate(value = value) %>% 
  ungroup()
params <- GAP_GEF_ratio %>% pull(param) %>% unique()
e.map <- e.map %>% 
  filter(mutant %in% GAP_GEF_ratio_mutants) %>% 
  filter(library_ORF %in% library_genes) %>% 
  #mutate("mutant" = factor(mutant, GAP_GEF_ratio_mutants)) %>% 
  arrange(mutant, library_ORF) %>% 
  mutate("score" = score)


rank_corr <- function(library_gene, param) {
  e.map.temp <- e.map[e.map$library_ORF == library_gene,]
  GAP_GEF_ratio_temp <- GAP_GEF_ratio[GAP_GEF_ratio$param == param, ]
  temp <- inner_join(e.map.temp, GAP_GEF_ratio_temp, by = "mutant")
  corr <- cor(temp$score, temp$value,
          use = "pairwise.complete.obs", method = "kendall")
}

pairs <- expand.grid(library_genes, params)
pairs <- list("library_gene" = as.character(pairs[,1]), "param" = as.character(pairs[,2]))
per_lib_gene_param_rank_correlation <- tibble(
  "library_ORF" = pairs$library_gene, 
  "parameter" = pairs$param,
  "kendall_rank_correlation" = vector("double", length(pairs$library_gene))
)
per_lib_gene_param_rank_correlation$kendall_rank_correlation <- pairs %>% pmap(rank_corr) %>% unlist()

per_lib_gene_param_rank_correlation <- per_lib_gene_param_rank_correlation %>% 
  inner_join(., orf_gene_index, by = c("library_ORF" = "ORF"))
per_lib_gene_param_rank_correlation %>% 
  ggplot(aes(x = kendall_rank_correlation)) + 
  geom_density() + 
  facet_wrap(~parameter)

gap_gef <- per_lib_gene_param_rank_correlation %>% 
  select(parameter, kendall_rank_correlation, gene_name) %>% 
  spread(parameter, kendall_rank_correlation) 

####### GAP kcat/Km and GEF kcat/Km
fit <- odregress(gap_gef$GAP_kcat_Km, gap_gef$GEF_kcat_Km)
cor <- signif(cor(gap_gef$GEF_kcat_Km, gap_gef$GAP_kcat_Km), digits = 2)
gap_gef %>% 
  ggplot(aes(x = GAP_kcat_Km, y = GEF_kcat_Km)) + #geom_point(alpha = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  ylab("relative GEF kcat/Km and E-MAP score") +
  xlab("relative GAP kcat/Km and E-MAP score") +
  ggtitle(str_c("Kendall rank correlation of Gsp1 mutants based on GAP or GEF relative efficiency and\n
  on genetic interactions (E-MAP S-score) with 1209 yeast genes - pearson correlation = ", cor)) +
  geom_abline(intercept = fit$coeff[2], slope = fit$coeff[1], size = 0.3) +
  geom_hline(yintercept = 0, size = 0.3) + geom_vline(xintercept = 0, size = 0.3)
ggsave("integrating_biophysical_parameters/GAP_vs_GEF_rank_cor_2ddensity_plot.pdf", width = 10, height = 8)

gap_gef %>% 
  ggplot(aes(x = GAP_kcat_Km, y = GEF_kcat_Km)) + geom_point(alpha = 0.5) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  ylab("relative GEF kcat/Km and E-MAP score") +
  xlab("relative GAP kcat/Km and E-MAP score") +
  ggtitle(str_c("Kendall rank correlation of Gsp1 mutants based on GAP or GEF relative efficiency and\n
                on genetic interactions (E-MAP S-score) with 1209 yeast genes - pearson correlation = ", cor)) +
  geom_abline(intercept = fit$coeff[2], slope = fit$coeff[1], size = 0.3) +
  geom_hline(yintercept = 0, size = 0.3) + geom_vline(xintercept = 0, size = 0.3)
ggsave("integrating_biophysical_parameters/GAP_vs_GEF_rank_cor_scatterplot.pdf", width = 10, height = 8)

# GAP kcat and GEF kcat
fit <- odregress(gap_gef$GAP_kcat, gap_gef$GEF_kcat)
cor <- signif(cor(gap_gef$GEF_kcat, gap_gef$GAP_kcat), digits = 2)
gap_gef %>% 
  ggplot(aes(x = GAP_kcat, y = GEF_kcat)) + #geom_point() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  ylab("GEF kcat and E-MAP score") +
  xlab("GAP kcat and E-MAP score") +
  ggtitle(str_c("Kendall rank correlation of Gsp1 mutants based on GAP and GEF kcat and\n
  on genetic interactions with 1536 yeast genes - cor = ", cor)) +
  geom_abline(intercept = fit$coeff[2], slope = fit$coeff[1], size = 0.3) +
  geom_hline(yintercept = 0, size = 0.3) + geom_vline(xintercept = 0, size = 0.3)
ggsave("integrating_biophysical_parameters/GAP_vs_GEF_kcat_rank_cor_2ddensity_plot.pdf", width = 10, height = 8)


fit <- odregress(gap_gef$GAP_Km, gap_gef$GEF_Km)
cor <- signif(cor(gap_gef$GEF_Km, gap_gef$GAP_Km), digits = 2)
gap_gef %>% 
  ggplot(aes(x = GAP_Km, y = GEF_Km)) + #geom_point() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  ylab("GEF kcat and E-MAP score") +
  xlab("GAP kcat and E-MAP score") +
  ggtitle(str_c("Kendall rank correlation of Gsp1 mutants based on GAP and GEF kcat and\n
                on genetic interactions with 1536 yeast genes - cor = ", cor)) +
  geom_abline(intercept = fit$coeff[2], slope = fit$coeff[1], size = 0.3) +
  geom_hline(yintercept = 0, size = 0.3) + geom_vline(xintercept = 0, size = 0.3)
ggsave("integrating_biophysical_parameters/GAP_vs_GEF_Km_rank_cor_2ddensity_plot.pdf", width = 10, height = 8)

test_complexes <- read_tsv("integrating_biophysical_parameters/test_complexes.txt", col_names = F) %>% 
  select("gene_name" = X2, "complex" = X3)

gap_gef_complex_eg <- gap_gef %>% 
  inner_join(., test_complexes, by = "gene_name")

#### kcat/Km
gap_gef %>% 
  ggplot(aes(x = GAP_kcat_Km, y = GEF_kcat_Km)) + 
  geom_point(alpha = 0.3) +
  geom_point(data = gap_gef_complex_eg, aes(color = complex), size = 4) +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  ylab("GEF kcat/Km and E-MAP score") +
  xlab("GAP kcat/Km and E-MAP score") +
  ggtitle("Kendall rank correlation of GTPase parameters and genetic interactions between yeast genes and Gsp1 mutants") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #scale_color_viridis(discrete = T)
  scale_color_manual( values = colors) + theme_gray()
ggsave("integrating_biophysical_parameters/GAP_vs_GEF_rel_eff_rank_cor_select_complexes_plot.pdf", height = 8, width = 10)


# #### kcat
# gap_gef %>% 
#   ggplot(aes(x = GAP_kcat, y = GEF_kcat)) + 
#   geom_point(alpha = 0.3) +
#   geom_point(data = gap_gef_complex_eg, aes(color = complex), size = 4) +
#   xlim(c(-1, 1)) + ylim(c(-1, 1)) +
#   ylab("GEF kcat and E-MAP score") +
#   xlab("GAP kcat and E-MAP score") +
#   ggtitle("Kendall rank correlation of GTPase parameters and genetic interactions between yeast genes and Gsp1 mutants") +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   scale_color_manual( values = colors) + theme_gray()
# ggsave("integrating_biophysical_parameters/GAP_vs_GEF_rel_kcat_rank_cor_select_complexes_plot.pdf", width = 14, height = 6)
# 
# 
# ### Km
# gap_gef %>% 
#   ggplot(aes(x = GAP_Km, y = GEF_Km)) + 
#   geom_point(alpha = 0.3) +
#   geom_point(data = gap_gef_complex_eg, aes(color = complex), size = 4) +
#   xlim(c(-1, 1)) + ylim(c(-1, 1)) +
#   ylab("GEF Km and E-MAP score") +
#   xlab("GAP Km and E-MAP score") +
#   ggtitle("Kendall rank correlation of GTPase parameters and genetic interactions between yeast genes and Gsp1 mutants") +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   scale_color_manual( values = colors) + theme_gray()
# ggsave("integrating_biophysical_parameters/GAP_vs_GEF_rel_Km_rank_cor_select_complexes_plot.pdf", width = 14, height = 6)
#cor(gap_gef$GEF_Km, gap_gef$GAP_Km)


# ### gamma2 and GAP_kcat_Km
# gap_gef %>% 
#   ggplot(aes(x = GAP_kcat_Km, y = gamma2)) + 
#   geom_point(alpha = 0.3) +
#   geom_point(data = gap_gef_complex_eg, aes(color = complex), size = 4) +
#   xlim(c(-1, 1)) + ylim(c(-1, 1)) +
#   ylab("% state2 and E-MAP score") +
#   xlab("GAP kcat/Km and E-MAP score") +
#   ggtitle("Kendall rank correlation of GTPase parameters and\n
#           genetic interactions between yeast genes and Gsp1 mutants") +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   scale_color_viridis(discrete = T)
# 
# cor(gap_gef$GAP_kcat_Km, gap_gef$gamma2, use = "pairwise.complete.obs")

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

ordered_genes <- gap_gef %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  gather(param, value, -gene_name) %>% 
  spread(gene_name, value) %>% 
  get_order_by_corr()

gap_gef %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  mutate("gene_name" = factor(gene_name, ordered_genes)) %>% 
  rename("GAP" = GAP_kcat_Km, "GEF" = GEF_kcat_Km) %>% 
  gather(param, value, -gene_name) %>% 
  ggplot(aes(param, gene_name, fill = value)) + geom_tile() +
  scale_fill_viridis() + 
  xlab("Kendall correlation between genetic interactions and GTPase cycle") +
  ylab("S. cerevisae gene") + theme(axis.text.y = element_blank())

write_tsv(gap_gef, "integrating_biophysical_parameters/kendall_correlation.txt")

##### ranking of library genes by correlations with either GAP or GEF relative efficiency
ranked_by_GAP <- gap_gef %>% 
  select(gene_name, GAP_kcat_Km) %>% 
  arrange(desc(GAP_kcat_Km))
write_tsv(ranked_by_GAP, "integrating_biophysical_parameters/ranked_by_GAP.rnk")
ranked_by_GEF <- gap_gef %>% 
  select(gene_name, GEF_kcat_Km) %>% 
  arrange(desc(GEF_kcat_Km))
write_tsv(ranked_by_GEF, "integrating_biophysical_parameters/ranked_by_GEF.rnk")


#### ranking on the diagonals
### shared correlations between GAP and GEF (ranking of genes that have correlations in the same direction)
### genes that are on the slope = 1 diagonal
shared_diagonal <- gap_gef %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  filter(sign(GAP_kcat_Km) == sign(GEF_kcat_Km)) %>% 
  mutate("projected_distance" = GAP_kcat_Km + GEF_kcat_Km) %>% 
  arrange(desc(projected_distance))
shared_kendall_corr_ranking <- shared_diagonal %>% select(gene_name, projected_distance)
write_tsv(shared_kendall_corr_ranking, "integrating_biophysical_parameters/GAP_GEF_shared_kendall_corr_ranking.rnk")
### genes that have orthogonal correlations between GAP and GEF (positive with one and negative with the other)
### genes that are on the slope = -1 diagonal
orthogonal_diagonal <- gap_gef %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  filter(sign(GAP_kcat_Km) != sign(GEF_kcat_Km)) %>% 
  mutate("projected_distance" = GAP_kcat_Km - GEF_kcat_Km) %>% 
  arrange(desc(projected_distance))
orthogonal_kendall_corr_ranking <- orthogonal_diagonal %>% select(gene_name, projected_distance)
write_tsv(orthogonal_kendall_corr_ranking, "integrating_biophysical_parameters/GAP_GEF_orthogonal_kendall_corr_ranking.rnk")

