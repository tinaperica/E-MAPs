# correlation of correlations, only APMS and filtered queries

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(psych)

load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')
name2ORF <- read_delim('spitzemap_name2ORF_index.txt', delim='\t')
apms_ORFs <- scan('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/SAINT_MS_hits.txt', what = 'character')
apms_ORFs <- c(apms_ORFs, 'YER009W') # add NTF2, which was not significantly pulled down in APMS, even by WT
apms_names <- filter(name2ORF, ORF %in% apms_ORFs) %>% pull(name)
GSP1_names <- filter(name2ORF, grepl('GSP1', name)) %>% pull(name)
core_partners <- c('srm1-g282s','srm1-ts','rna1-s116f','rna1-1','srp1-5001',
                   'yrb1-51','kap95-e126k','mog1','ntf2-5001','ntf2-h104y')


# start <- Sys.time()

# correlation matrix for all significant ORFs, including Gsp1 (takes 3.24 min)
# cormat <-
#   filtered_correlations %>%
#   filter(two_sided_p_value < 0.05) %>%
#   filter(!grepl('GSP1', query_uniq2)) %>% 
#   select(query_uniq1, query_uniq2, pearson) %>%
#   spread(query_uniq2, pearson) %>%
#   column_to_rownames('query_uniq1') %>%
#   as.matrix() %>%
#   t() %>% 
#   # cor(., use = "pairwise.complete.obs", method = "pearson")
#   corr.test(., use = "pairwise", method = "pearson", ci=F)
# 
# Sys.time() - start
# 
# save(cormat, file = 'corr_of_corr.RData')

load('corr_of_corr.RData')

#cormat$p has the adjusted pvalues in the upper tri
adj_pmask <- cormat$p
adj_pmask[lower.tri(adj_pmask)] = 0
adj_pmask = adj_pmask + t(adj_pmask)
adj_pmask <- 1*(adj_pmask < 0.05)
cormat_filtered_corrected <- as.matrix(cormat$r)*adj_pmask

#cormat$p has the unadjusted pvalues in the lower tri
pmask <- cormat$p
pmask[upper.tri(pmask)] = 0
pmask = pmask + t(pmask)
pmask <- 1*(pmask < 0.05)
cormat_filtered <- as.matrix(cormat$r)*pmask


GSP1_to_keep <- intersect(rownames(cormat_filtered), GSP1_names)
columns_to_keep <- intersect(colnames(cormat_filtered), apms_names)
columns_to_keep_noGsp1 <- setdiff(intersect(colnames(cormat_filtered), apms_names), GSP1_to_keep)

pdf('all_vs_gsp1_filtered_uncorrected.pdf', height = 7, width = 25)
Heatmap(cormat_filtered[GSP1_to_keep, columns_to_keep_noGsp1],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('all_vs_all_filtered_uncorrected.pdf', height = 25, width = 25)
Heatmap(cormat_filtered[columns_to_keep, columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('all_vs_gsp1_filtered_corrected.pdf', height = 7, width = 25)
Heatmap(cormat_filtered_corrected[GSP1_to_keep, columns_to_keep_noGsp1],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('all_vs_all_filtered_corrected.pdf', height = 25, width = 25)
Heatmap(cormat_filtered_corrected[columns_to_keep, columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('all_vs_gsp1_unfiltered.pdf', height = 7, width = 25)
Heatmap(cormat$r[GSP1_to_keep, columns_to_keep_noGsp1],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('all_vs_all_unfiltered.pdf', height = 25, width = 25)
Heatmap(cormat$r[columns_to_keep, columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()


# Tried to run same, but with *only% Gsp1 mutants as correlation vectors
# But there aren't enough complete cases for this to work.
cormat2 <-
  filtered_correlations %>%
  filter(two_sided_p_value < 0.05) %>%
  filter(grepl('GSP1', query_uniq2)) %>% 
  select(query_uniq1, query_uniq2, pearson) %>%
  spread(query_uniq2, pearson) %>%
  column_to_rownames('query_uniq1') %>%
  as.matrix() %>%
  t() %>% 
  # cor(., use = "pairwise.complete.obs", method = "pearson")
  corr.test(., use = "pairwise", method = "pearson", ci=F)

columns_to_keep <- intersect(colnames(cormat2$r), apms_names)

pdf('all_vs_all_onlyGsp1_columns.pdf', height = 25, width = 25)
Heatmap(cormat2$r[columns_to_keep, columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

