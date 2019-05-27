# correlation of correlations, only APMS and filtered queries

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')
name2ORF <- read_delim('spitzemap_name2ORF_index.txt', delim='\t')
apms_ORFs <- scan('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/SAINT_MS_hits.txt', what = 'character')
apms_ORFs <- c(apms_ORFs, 'YER009W') # add NTF2, which was not significantly pulled down in APMS, even by WT
apms_names <- filter(name2ORF, ORF %in% apms_ORFs) %>% pull(name)
GSP1_names <- filter(name2ORF, grepl('GSP1', name)) %>% pull(name)
core_partners <- c('srm1-g282s','srm1-ts','rna1-s116f','rna1-1','srp1-5001',
                   'yrb1-51','kap95-e126k','mog1','ntf2-5001','ntf2-h104y')


# correlation matrix for just apms_ORFs, including Gsp1
cormat <- 
  filtered_correlations %>% 
  left_join(name2ORF, by = c('query_uniq1' = 'name')) %>% 
  rename('ORF1' = ORF) %>% 
  left_join(name2ORF, by = c('query_uniq2' = 'name')) %>% 
  rename('ORF2' = ORF) %>% 
  filter(ORF1 %in% apms_ORFs, ORF2 %in% apms_ORFs) %>% 
  filter(!grepl('GSP1',query_uniq2)) %>% 
  filter(two_sided_p_value < 0.05) %>% 
  select(query_uniq1, query_uniq2, pearson) %>%
  spread(query_uniq2, pearson) %>%
  column_to_rownames('query_uniq1') %>% 
  as.matrix() %>% 
  cor(., use = "pairwise.complete.obs", method = "pearson")

# Full Heatmap
pdf('corr_of_corr_APMSonly_withGsp1_sym.pdf', height = 25, width = 25)
Heatmap(cormat, name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# Heatmap of just Gsp1 vs APMS partners
pdf('corr_of_corr_APMSonly_withGsp1_assym.pdf', height = 6, width = 25)
cormat %>% 
  data.frame() %>% 
  mutate('query_uniq1' = rownames(cormat)) %>% 
  select(query_uniq1, everything()) %>% 
  filter(grepl('GSP1', query_uniq1)) %>% 
  select(-contains('GSP1')) %>% 
  column_to_rownames('query_uniq1') %>% 
  as.matrix() %>% 
  Heatmap(name = 'Corr. of corr.',
          col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
          heatmap_legend_param = 'bottom',
          row_title = 'query', column_title = 'query',
          clustering_method_rows = 'complete', clustering_method_columns = 'complete',
          row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
          row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# correlation matrix for just apms_ORFs, excluding Gsp1
cormat_woutGsp1 <-
  filtered_correlations %>% 
  left_join(name2ORF, by = c('query_uniq1' = 'name')) %>% 
  rename('ORF1' = ORF) %>% 
  left_join(name2ORF, by = c('query_uniq2' = 'name')) %>% 
  rename('ORF2' = ORF) %>% 
  filter(ORF1 %in% apms_ORFs, ORF2 %in% apms_ORFs) %>% 
  filter(two_sided_p_value < 0.05) %>% 
  filter(!grepl('GSP1', query_uniq1), !grepl('GSP1', query_uniq2)) %>%  
  select(query_uniq1, query_uniq2, pearson) %>%
  spread(query_uniq2, pearson) %>%
  column_to_rownames('query_uniq1') %>% 
  as.matrix() %>% 
  cor(., use = "pairwise.complete.obs", method = "pearson")

pdf('corr_of_corr_APMSonly_woutGsp1.pdf', height = 25, width = 25)
Heatmap(cormat_woutGsp1,
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# correlation matrix for all significant ORFs, including Gsp1 (takes 3.24 min)
cormat_all_sig_queries <-
  filtered_correlations %>%
  filter(two_sided_p_value < 0.05) %>%
  filter(!grepl('GSP1',query_uniq2)) %>% 
  select(query_uniq1, query_uniq2, pearson) %>%
  spread(query_uniq2, pearson) %>%
  column_to_rownames('query_uniq1') %>%
  as.matrix() %>%
  cor(., use = "pairwise.complete.obs", method = "pearson")

columns_to_keep <- intersect(colnames(cormat_all_sig_queries), apms_names)
GSP1_to_keep <- intersect(colnames(cormat_all_sig_queries), GSP1_names)

# plot the correlation matrix of just APMS partners and Gsp1, but corr of corrs computed with all queries
pdf('corr_of_corr_allQs_withGsp1_sym.pdf', height = 25, width = 25)
Heatmap(cormat_all_sig_queries[columns_to_keep,columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# plot the correlation matrix of just APMS partners and Gsp1, but corr of corrs computed with all queries
pdf('corr_of_corr_allQs_withGsp1_assym.pdf', height = 6, width = 25)
Heatmap(cormat_all_sig_queries[GSP1_to_keep,columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# correlation matrix for all significant ORFs, excluding Gsp1 (takes 3.24 min)
# cormat_all_sig_queries_woutGsp1 <- 
#   filtered_correlations %>% 
#   filter(two_sided_p_value < 0.05) %>% 
#   filter(!grepl('GSP1', query_uniq1), !grepl('GSP1', query_uniq2)) %>%  
#   select(query_uniq1, query_uniq2, pearson) %>%
#   spread(query_uniq2, pearson) %>%
#   column_to_rownames('query_uniq1') %>% 
#   as.matrix() %>% 
#   cor(., use = "pairwise.complete.obs", method = "pearson")

columns_to_keep <- intersect(colnames(cormat_all_sig_queries_woutGsp1), apms_names)

# plot the correlation matrix of just APMS partners, excluding Gsp1
pdf('corr_of_corr_allQs_woutGsp1_sym.pdf', height = 25, width = 25)
Heatmap(cormat_all_sig_queries_woutGsp1[columns_to_keep,columns_to_keep],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()




# plot small correlation matrices of just core partners
pdf('partners_small_corr_of_corr_heatmaps/corr_of_corr_APMSonly_withGsp1_partners.pdf', height = 5, width = 5)
Heatmap(cormat[core_partners,core_partners],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(10, "mm"), column_dend_height = unit(10, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('partners_small_corr_of_corr_heatmaps/corr_of_corr_APMSonly_woutGsp1_partners.pdf', height = 5, width = 5)
Heatmap(cormat_woutGsp1[core_partners,core_partners],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(10, "mm"), column_dend_height = unit(10, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('partners_small_corr_of_corr_heatmaps/corr_of_corr_allQs_withGsp1_partners.pdf', height = 5, width = 5)
Heatmap(cormat_all_sig_queries[core_partners,core_partners],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(10, "mm"), column_dend_height = unit(10, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

pdf('partners_small_corr_of_corr_heatmaps/corr_of_corr_allQs_woutGsp1_partners.pdf', height = 5, width = 5)
Heatmap(cormat_all_sig_queries_woutGsp1[core_partners,core_partners],
        name = 'Corr. of corr.',
        col = colorRamp2(c(-1, 0, 1), c("orange", "white", "purple")),
        heatmap_legend_param = 'bottom',
        row_title = 'query', column_title = 'query',
        clustering_method_rows = 'complete', clustering_method_columns = 'complete',
        row_dend_width = unit(10, "mm"), column_dend_height = unit(10, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

