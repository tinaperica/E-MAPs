
library(tidyverse)
library(gplots)
library(ComplexHeatmap)
library(circlize)

# define clustering function
clustfn <- function(mat) {
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  
  n_NA = length(which(sapply(dissim, function(x)all(is.na(x)))))
  if ((n_NA) > 0) {
    print(paste0('setting ', n_NA, ' dissimilarity value(s) to 1'))
    dissim[which(sapply(dissim, function(x)all(is.na(x))))] <- 1.0 # if dissim is NA, set to max distance
  }
  hc <- hclust(dissim, method = "average")
  return(hc)
}

# cluster full emap
mat <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
row.hc <- clustfn(t(mat))
col.hc <- clustfn(mat)
full_emap_arrays.hc <- col.hc # save for later printing

# make a heatmap of the full emap, with mutants and arrays clustered
pdf('images/full_emap.pdf', height = 10, width = 60)
Heatmap(mat, name = 'S-score',
        col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
        heatmap_legend_param = 'bottom',
        row_title = 'gsp1 mutant',
        column_title = 'Full EMAP, clustered by pearson correlation',
        column_title_gp = gpar(fontsize = 28),
        cluster_rows = row.hc, cluster_columns = col.hc,
        row_dend_reorder = FALSE, column_dend_reorder = FALSE,
        row_dend_width = unit(50, "mm"), column_dend_height = unit(50, "mm"),
        column_names_gp = gpar(fontsize = 3))
dev.off()

# also plot emap without mutant clustering (previously ordered by residue)
pdf('images/full_emap_by_res.pdf', height = 10, width = 60)
Heatmap(mat, name = 'S-score',
        col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
        row_title = 'gsp1 mutant',
        column_title = 'Full EMAP, mutants ordered by residue',
        column_title_gp = gpar(fontsize = 28),
        cluster_rows = FALSE, cluster_columns = col.hc,
        row_dend_reorder = FALSE, column_dend_reorder = FALSE,
        row_dend_width = unit(50, "mm"), column_dend_height = unit(50, "mm"),
        column_names_gp = gpar(fontsize = 3))
dev.off()

# get list of mutants for which we have GAP and GEF data
mutants_with_GAP_GEF_data_ordered <-
  read_delim('datasets/GAP_GEF_kinetics_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  arrange(GAP_GEF_ratio) %>%
  pull(mutant)

# subset emap to only include mutants with GAP GEF data
mat <- mat[mutants_with_GAP_GEF_data_ordered,]
row.hc <- clustfn(t(mat))
col.hc <- clustfn(mat)

# make heatmap
pdf('images/full_emap_GAPGEF_only.pdf', height = 10, width = 60)
Heatmap(mat, name = 'S-score',
        col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
        row_title = 'gsp1 mutant',
        column_title = 'Full EMAP, GAP GEF complete clustered by pearson correlation',
        column_title_gp = gpar(fontsize = 28),
        cluster_rows = row.hc, cluster_columns = col.hc,
        row_dend_reorder = FALSE, column_dend_reorder = FALSE,
        row_dend_width = unit(50, "mm"), column_dend_height = unit(50, "mm"),
        column_names_gp = gpar(fontsize = 3))
dev.off()

# Read in interface ΔrASA for clustering
mat <- read_delim('datasets/core_deltarASA_by_mutant_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
row.hc <- clustfn(t(mat))
col.hc <- clustfn(mat)

# heatmap
pdf('images/interface_clustering.pdf', height = 10, width = 8)
Heatmap(mat, name = 'Delta rASA upon binding',
        col = colorRamp2(c(0,1), c("black", "cyan")),
        row_title = 'gsp1 mutant',
        column_title = 'Mutants clustered by interface',
        column_title_gp = gpar(fontsize = 18),
        cluster_rows = row.hc, cluster_columns = col.hc,
        row_dend_reorder = FALSE, column_dend_reorder = FALSE,
        row_dend_width = unit(20, "mm"), column_dend_height = unit(20, "mm"))
dev.off()

# make a heatmap of the full emap, but mutants ordered by interface

mat <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  filter(! mutant %in% c('GSP1-NAT', 'NTER3XFLAG WT', 'CTER3XFLAG WT')) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
col.hc <- clustfn(mat)

pdf('images/full_emap_by_interface.pdf', height = 10, width = 60)
Heatmap(mat, name = 'S-score',
        col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
        heatmap_legend_param = 'bottom',
        row_title = 'gsp1 mutant',
        column_title = 'Full EMAP, mutants clustered by interface',
        column_title_gp = gpar(fontsize = 28),
        cluster_rows = row.hc, cluster_columns = full_emap_arrays.hc,
        row_dend_reorder = FALSE, column_dend_reorder = FALSE,
        row_dend_width = unit(50, "mm"), column_dend_height = unit(50, "mm"),
        column_names_gp = gpar(fontsize = 3))
dev.off()


# order GAP/GEF dataset
mat <- read_delim('datasets/GAP_GEF_kinetics_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  arrange(GAP_GEF_ratio) %>%
  mutate('delta' = ln_GAP_kcat_Km - ln_GEF_kcat_Km) %>% 
  column_to_rownames(var="mutant") %>%
  select(ln_GAP_kcat_Km, ln_GEF_kcat_Km, delta) %>%
  as.matrix()

# make a heatmap of the GAP GEF values, ordered by ratio of relative efficiency (GAP/GEF)
pdf('images/GAPGEF_clustering.pdf', height = 10, width = 6)
Heatmap(mat, name = 'ln(kcat/Km)',
        col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
        row_title = 'gsp1 mutant',
        column_title = 'Mutants ordered by ratio of GAP, GEF efficiency',
        column_title_gp = gpar(fontsize = 12),
        cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()
