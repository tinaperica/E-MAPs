# correlation of correlations, only APMS and filtered queries

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(psych)

load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')
name2ORF <- read_delim('cjm_corr/spitzemap_name2ORF_index.txt', delim='\t')
apms_ORFs <- scan('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/SAINT_MS_hits.txt', what = 'character')
apms_ORFs <- c(apms_ORFs, 'YER009W') # add NTF2, which was not significantly pulled down in APMS, even by WT
apms_names <- filter(name2ORF, ORF %in% apms_ORFs) %>% pull(name)
GSP1_names <- filter(name2ORF, grepl('GSP1', name)) %>% pull(name)
core_partners <- c('srm1-g282s','srm1-ts','rna1-s116f','rna1-1','srp1-5001',
                   'yrb1-51','kap95-e126k','mog1','ntf2-5001','ntf2-h104y')


filter_out_na <- function(.mat) {
  na_count <- apply(.mat, 1, function(x) sum(is.na(x)))
  .mat[which(ncol(.mat) - na_count >= 3), which(ncol(.mat) - na_count >= 3)]
}

plot_corr_of_corr <- function(cormat, name, rm.NA = T) {
  
  cormat_unfiltered <- cormat$r
  cormat_unfiltered <- filter_out_na(cormat_unfiltered)
  
  #cormat$p has the unadjusted pvalues in the lower tri
  pmask <- cormat$p
  pmask[upper.tri(pmask)] = 0
  pmask = pmask + t(pmask)
  pmask <- 1*(pmask < 0.05)
  if(rm.NA == F) {
    pmask[pmask == 0] <- NA
  }
  cormat_filtered <- as.matrix(cormat$r)*pmask
  comrat_filtered <- filter_out_na(cormat_filtered)
  
  
  #cormat$p has the adjusted pvalues in the upper tri
  adj_pmask <- cormat$p
  adj_pmask[lower.tri(adj_pmask)] = 0
  adj_pmask = adj_pmask + t(adj_pmask)
  adj_pmask <- 1*(adj_pmask < 0.05)
  if(rm.NA == F) {
    adj_pmask[adj_pmask == 0] <- NA
  }
  cormat_filtered_corrected <- as.matrix(cormat$r)*adj_pmask
  cormat_filtered_corrected <- filter_out_na(cormat_filtered_corrected)
  
  GSP1_to_keep <- intersect(rownames(cormat_unfiltered), GSP1_names)
  columns_to_keep <- intersect(colnames(cormat_unfiltered), apms_names)
  columns_to_keep_noGsp1 <- setdiff(columns_to_keep, GSP1_to_keep)
  
  colramp = colorRamp2(c(-0.5, 0, 0.5), c("orange", "white", "purple"))
  pdf(name, height = 7, width = 25)

  # try(print(Heatmap(cormat_unfiltered[columns_to_keep_noGsp1, columns_to_keep_noGsp1],
  #                   col = colramp,
  #                   column_title = 'Correlations of correlations, unfiltered',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
  #                   row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))))
  # 
  # try(print(Heatmap(cormat_unfiltered[GSP1_to_keep, columns_to_keep_noGsp1], col = colramp,
  #                   column_title = 'Correlations of correlations, unfiltered',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(30, "mm"), column_dend_height = unit(150, "mm"),
  #                   row_names_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 8))))
  # 
  # try(print(Heatmap(cormat_unfiltered[GSP1_to_keep, core_partners], col = colramp,
  #                   column_title = 'Correlations of correlations, unfiltered',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(150, "mm"), column_dend_height = unit(150, "mm"),
  #                   row_names_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 30))))
  
  GSP1_to_keep <- intersect(rownames(cormat_filtered), GSP1_names)
  columns_to_keep <- intersect(colnames(cormat_filtered), apms_names)
  columns_to_keep_noGsp1 <- setdiff(columns_to_keep, GSP1_to_keep)
  # 
  # try(print(Heatmap(cormat_filtered[columns_to_keep_noGsp1, columns_to_keep_noGsp1],
  #               col = colramp,
  #               column_title = 'Correlations of correlations, p < 0.05',
  #               column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #               row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
  #               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))))
  
  try(print(Heatmap(cormat_unfiltered[GSP1_to_keep, union(core_partners, GSP1_names)], col = colramp,
                    column_title = 'Correlations of correlations, p < 0.05',
                    column_title_gp = gpar(fontsize = 40, fontface = "bold"),
                    row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
                    row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 8))))
  
  # try(print(Heatmap(cormat_filtered[GSP1_to_keep, core_partners], col = colramp,
  #                   column_title = 'Correlations of correlations, p < 0.05',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(150, "mm"), column_dend_height = unit(150, "mm"),
  #                   row_names_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 30))))
  
  # GSP1_to_keep <- intersect(rownames(cormat_filtered_corrected), GSP1_names)
  # columns_to_keep <- intersect(colnames(cormat_filtered_corrected), apms_names)
  # columns_to_keep_noGsp1 <- setdiff(columns_to_keep, GSP1_to_keep)
  # 
  # try(print(Heatmap(cormat_filtered_corrected[columns_to_keep_noGsp1, columns_to_keep_noGsp1], col = colramp,
  #                   column_title = 'Correlations of correlations, adjusted p < 0.05',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
  #                   row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))))
  # 
  # try(print(Heatmap(cormat_filtered_corrected[GSP1_to_keep, columns_to_keep_noGsp1], col = colramp,
  #                   column_title = 'Correlations of correlations, adjusted p < 0.05',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(30, "mm"), column_dend_height = unit(150, "mm"),
  #                   row_names_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 8))))
  # 
  # try(print(Heatmap(cormat_filtered_corrected[GSP1_to_keep, core_partners], col = colramp,
  #                   column_title = 'Correlations of correlations, adjusted p < 0.05',
  #                   column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  #                   row_dend_width = unit(150, "mm"), column_dend_height = unit(150, "mm"),
  #                   row_names_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 30))))
  
  dev.off()
}


pdf('hist_pval_corrs.pdf')
x <- filtered_correlations %>% 
  select(query_uniq1, query_uniq2, pearson, two_sided_p_value) %>% 
  mutate('sig' = case_when((two_sided_p_value > 0.05) ~ FALSE, TRUE ~ TRUE))
ggplot(x, aes(x=pearson)) +
  geom_histogram(data = filter(x, sig == F), fill = 'red', alpha = 0.2, bins = 5000) +
  geom_histogram(data = filter(x, sig == T), fill = 'blue', alpha = 0.2, bins = 5000) +
  xlim(c(-0.25, 0.25)) + ggtitle('p-values of correlations of S-scores')
dev.off()

pdf('srm1.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','srm1-ts')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `srm1-ts`)) + geom_point()
dev.off()

pdf('srm1_rna1.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','rna1-1')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `rna1-1`)) + geom_point()
dev.off()

pdf('srm1_acf4.pdf')
filtered_correlations %>% 
  select(starts_with('query_uniq'), pearson) %>% 
  filter(query_uniq2 %in% c('srm1-g282s','acf4')) %>%
  spread(query_uniq2, pearson) %>%
  ggplot(aes(x = `srm1-g282s`, y = `acf4`)) + geom_point()
dev.off()



run_fn <- function(name, filter_pval = T, use_Gsp1 = F, use_weighted = F, load=T) {
  
  start <- Sys.time()
  
  if (load) {
    # load('corr_of_weighted_corr_unfiltered.RData')
    load(paste0(name, '.RData'))
  } else {
    
    df <-
      filtered_correlations %>%
      filter(
        if (filter_pval) {
          two_sided_p_value < 0.05
        } else {
          two_sided_p_value == two_sided_p_value
        }) %>%
      filter(
        if (use_Gsp1) {
          grepl('GSP1', query_uniq2)
          } else {
          !grepl('GSP1', query_uniq2)
        })
      
    if(use_weighted) {
      cormat <-
        df %>% 
        select(query_uniq1, query_uniq2, min_weighted_pearson) %>% 
        spread(query_uniq2, min_weighted_pearson) %>%
        column_to_rownames('query_uniq1') %>% as.matrix() %>% t() %>%
        corr.test(., use = "pairwise", method = "pearson", ci=F)
    } else {
      cormat <-
        df %>% 
        select(query_uniq1, query_uniq2, pearson) %>% 
        spread(query_uniq2, pearson) %>%
        column_to_rownames('query_uniq1') %>% as.matrix() %>% t() %>%
        corr.test(., use = "pairwise", method = "pearson", ci=F)
    }
    save(cormat, file = paste0(name, '.RData'))
  }
  plot_corr_of_corr(cormat, paste0(name, '2.pdf'))
  
  print(Sys.time() - start)
  
}


run_fn('corr_of_unweighted_corr_filtered', use_weighted = F, filter_pval = T, use_Gsp1 = F)




run_fn('corr_of_unweighted_corr_unfiltered', use_weighted = F, filter_pval = F, use_Gsp1 = F)
run_fn('corr_of_unweighted_corr_filtered', use_weighted = F, filter_pval = T, use_Gsp1 = F)
run_fn('corr_of_unweighted_corr_unfiltered_Gsp1only', use_weighted = F, filter_pval = F, use_Gsp1 = T)
run_fn('corr_of_unweighted_corr_filtered_Gsp1only', use_weighted = F, filter_pval = T, use_Gsp1 = T)
run_fn('corr_of_weighted_corr_unfiltered', use_weighted = T, filter_pval = F, use_Gsp1 = F)
run_fn('corr_of_weighted_corr_filtered', use_weighted = T, filter_pval = T, use_Gsp1 = F)
run_fn('corr_of_weighted_corr_unfiltered_Gsp1only', use_weighted = T, filter_pval = F, use_Gsp1 = T)
run_fn('corr_of_weighted_corr_filtered_Gsp1only', use_weighted = T, filter_pval = T, use_Gsp1 = T)


