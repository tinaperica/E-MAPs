#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(cowplot)
library(viridis)
library(scales)

# GAP_GEF_complete_only = args[1] # include mutants without GAP/GEF kinetic data
# use_correlations <- args[2] # score or correlations
GAP_GEF_complete_only = FALSE
use_correlations = TRUE

# load emap
if (use_correlations) {
  emap <-
    read_delim('datasets/correlations_filtered_by_query_and_corr_strength.txt', delim = '\t', col_types = cols()) %>%
    separate('query_uniq1', into = c('gsp1','mutant'), sep = ' - ') %>%
    select(-gsp1) %>% 
    gather(name, score, -mutant) %>% 
    left_join(read_delim('datasets/spitzemap_name2ORF_index.txt', delim = '\t'), by = 'name') %>% 
    mutate("mutant" = case_when(mutant == "CTER3XFLAG WT" ~ "C_WT",
                                mutant == "NTER3XFLAG WT" ~ "N_WT",
                                TRUE ~ mutant))
} else {
  emap <-
    bind_cols(
      read_delim('datasets/gsp1_pEMAP_avg_merged_gene_names.txt', delim='\t', col_types=cols()) %>% 
        gather('name', 'score', -Gene),
      read_delim('datasets/gsp1_pEMAP_avg_merged.txt', delim='\t', col_types=cols()) %>% 
        gather('ORF', 'score', -Gene)) %>% 
    separate(Gene, into = c('gsp1', 'mutant'), sep = ' - ', extra = 'merge') %>% 
    select(mutant, name, ORF, score) %>% 
    mutate("mutant" = case_when(mutant == "GSP1-NAT" ~ "WT",
                                mutant == "CTER3XFLAG WT" ~ "C_WT",
                                mutant == "NTER3XFLAG WT" ~ "N_WT",
                                TRUE ~ mutant))
}

# GAP/GEF ratio
GAP_GEF_kinetics <- read_delim('datasets/GAP_GEF_kinetics_for_clustering.txt', delim = '\t', col_types = cols())
GAP_GEF_ratio_mutants <- GAP_GEF_kinetics %>% pull(mutant) %>% unique()
GAP_GEF_with_NA <-
  left_join(unique(select(emap, mutant)), GAP_GEF_kinetics, by = 'mutant')
if (GAP_GEF_complete_only) {
  GAP_GEF_with_NA <-
    GAP_GEF_with_NA %>%
    replace_na(list(GAP_kcat_Km = 0, GEF_kcat_Km = 0, GAP_GEF_ratio = 1,
                    ln_GAP_kcat_Km = 0, ln_GEF_kcat_Km = 0))
  emap <- filter(emap, (mutant %in% GAP_GEF_ratio_mutants) | (mutant == 'WT'))
}

# read in apms data for point plot
apms_data <-
  read_tsv('datasets/apms_log2_fold_change.txt', col_types = cols()) %>% 
  filter(gene_name %in% c("SRM1", "RNA1") & norm == "eqM") %>% 
  select(mutant, log2FC, adj.pvalue, sample, gene_name)
apms_data <-
  apms_data %>% 
  select(-mutant, -adj.pvalue) %>% 
  spread(gene_name, log2FC) %>% 
  mutate('diffGAPGEF' = RNA1 - SRM1) %>% 
  select(sample, diffGAPGEF) %>% 
  gather(key = gene_name, value = log2FC, diffGAPGEF) %>% 
  bind_rows(apms_data) %>% 
  mutate(adj.pvalue = ifelse(is.na(adj.pvalue), 0, adj.pvalue)) %>% 
  separate(sample, into = c('tag','mutant'), remove=T) %>% 
  unite(partner_tag, tag, gene_name)
missing_apms_mutants <-
  setdiff(unique(select(emap, mutant)),
          unique(select(apms_data, mutant))) %>% pull()
missing_mutants_df <-
  expand.grid(missing_apms_mutants,
              unique(apms_data$partner_tag))
colnames(missing_mutants_df) <- c('mutant','partner_tag')
apms_with_NA <- bind_rows(apms_data, missing_mutants_df)

# Read in interface ΔrASA for clustering
interface <-
  emap %>% 
  select(mutant) %>% 
  mutate(yeastresnum = as.numeric(substr(mutant, 2, nchar(mutant)-1))) %>% 
  left_join(read_delim('../../ran_structures/SASA_interfaces.txt',
                       delim = '\t', col_types = cols()),
            by = 'yeastresnum') %>% 
  select(mutant, partner, deltarASA) %>%
  unique() %>%
  spread(partner, deltarASA, fill=FALSE) %>%
  select(-`<NA>`) %>%
  gather(partner, deltarASA, -mutant)

## if using scores (not correlations), we want to keep only array genes
# that have more than one member that has a significant E-MAP score
# with at least one mutant
if (use_correlations) {
  ORFs_to_keep <-
    emap %>% 
    pull(ORF) %>% 
    unique()
} else {
  ORFs_to_keep <-
    emap %>% 
    filter(score < -2 | score > 2) %>%
    pull(ORF) %>% 
    unique()
}

emap_data <-
  read_delim('gene_sets/gene_sets_combined.txt', delim='\t', col_types = cols()) %>% 
  filter(ORF %in% ORFs_to_keep) %>% 
  group_by(gene_set) %>% 
  mutate(n_members = n()) %>% 
  ungroup() %>% 
  filter(n_members > 1) %>% 
  inner_join(emap, by = 'ORF') %>% 
  select(-ORF, -n_members, -Gene)

gene_set_list <-
  emap_data %>% 
  select(gene_set, set_label) %>% 
  unique()

plot_list <- list()

for (i in seq(nrow(gene_set_list))) {

  gset <- gene_set_list$gene_set[i]
  gset_label <- gene_set_list$set_label[i]
  print(paste(i, gset, gset_label, sep='--'))
  
  temp.emap <-
    emap_data %>% 
    filter(gene_set == gset, set_label == gset_label) %>% 
    select(-gene_set, -set_label) %>%
    unique() %>% 
    spread(name, score)
  
  # make matrix for distance calculation
  mat <-
    temp.emap %>% 
    column_to_rownames('mutant') %>% 
    as.matrix()
  mat[is.na(mat)] <- 0
  distance <- dist(mat, method = "euclidean")
  
  # use multi-dimensional scaling to order the mutants by their first PC
  ordered_genes <-
    cmdscale(distance, eig = T, k = 1) %>% 
    `$`(point) %>% 
    as_tibble(rownames = 'mutant') %>% 
    arrange(V1) %>% 
    pull(mutant)
  
  # heatmap of emap S-scores
  if (use_correlations) {
    scale_limit <- c(-0.5,0.5)
    low_color <- 'red'
    mid_color <- 'white'
    high_color <- 'green'
    scatter_label <- 'correlation'
  } else {
    scale_limit <- c(-3,3)
    low_color <- 'blue'
    mid_color <- 'black'
    high_color <- 'yellow'
    scatter_label <- 'EMAP score'
  }
  
  if (ncol(temp.emap > 12)) {
    label_size = 6
  } else {
    label_size = 10
  }
  
  emap_hm <-
    temp.emap %>% 
    gather(gene, score, -mutant) %>% 
    mutate('mutant' = factor(mutant, ordered_genes)) %>% 
    ggplot(aes(mutant, gene)) +
    geom_tile(aes(fill = score), color = 'white') +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color,
                         midpoint = 0, limit=scale_limit, oob=squish,
                         labels = trans_format('identity', function(x) round(x, 2))) +
    theme_grey() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = label_size),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(0.9, 1.04),
          legend.direction = 'horizontal',
          legend.key.size = unit(0.5, "cm")) +
    ggtitle(paste(gset, gset_label, sep='\n'))
  
  # compute mean WT emap S-score to place on the scatter plot
  mean_wt_score <- emap_data %>% 
    filter(gene_set == gset & mutant == "WT") %>%
    summarize('mean' = mean(score)) %>% pull(mean)
  
  # scatter plot of emap S-scores
  emap_plot <- temp.emap %>% 
    gather(key = gene, value = score, -mutant) %>% 
    mutate('mutant' = factor(mutant, ordered_genes)) %>% 
    ggplot(aes(x = mutant, y = score, color = gene)) + 
    geom_point(size = 4, alpha = 0.75) +
    geom_hline(yintercept = mean_wt_score) +
    theme_grey() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = c(0.5, 0.2),
          legend.direction = 'horizontal',
          legend.key.size = unit(0.2, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_viridis(discrete = T) +
    guides(colour = guide_legend(nrow = ceiling(ncol(temp.emap)/8))) +
    labs(y = scatter_label)
  
  # prepare GAP, GEF kinetics data
  GAP_GEF_for_plotting <- mutate(GAP_GEF_with_NA, 'mutant' = factor(mutant, ordered_genes))
  if(GAP_GEF_complete_only) {
    GAP_GEF_for_plotting <- drop_na(GAP_GEF_for_plotting)
  }
  
  # bar plot of GAP kcat/Km
  GAP_plot <-
    ggplot(GAP_GEF_for_plotting, aes(x = mutant, y = GAP_kcat_Km)) +
    geom_bar(stat='identity') +
    theme_grey() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(y = 'GAP kcat/Km\nMUT/WT')

  # bar plot of GEF kcat/Km
  GEF_plot <-
    ggplot(GAP_GEF_for_plotting, aes(x = mutant, y = GEF_kcat_Km)) +
    geom_bar(stat='identity') +
    theme_grey() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(y = 'GEF kcat/Km\nMUT/WT')
  
  # bar plot of GAP, GEF kcat/Km ratio
  GAPGEF_plot <-
    ggplot(GAP_GEF_for_plotting, aes(x = mutant, y = log(GAP_GEF_ratio))) +
    geom_bar(stat='identity') +
    theme_grey() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(y = 'ln(GAP/GEF)') + geom_hline(yintercept=0) 
  
  # point plot of APMS GAP GEF log2FC, including GAP-GEF difference
  apms_plot <-
    apms_with_NA %>% 
    mutate("mutant" = factor(mutant, ordered_genes)) %>%
    arrange(mutant) %>% 
    ggplot(aes(x = mutant, y = partner_tag, fill = log2FC, size = adj.pvalue)) +
    geom_point(shape = 21, stroke = 0.1) +
    scale_fill_gradient2() + 
    scale_size("adj.pvalue", range = c(4, 0.1),
               breaks = c(0, 0.001, 0.0375, 0.05, 0.1, 0.25, 0.5)) +
    ylab('APMS log2FC') +
    theme_grey() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "top", legend.direction = 'horizontal',
          legend.key.size = unit(0.3, "cm"))

  # heatmap of interface deltarASA
  iface_plot <-
    interface %>% 
    mutate('mutant' = factor(mutant, ordered_genes)) %>% 
    ggplot(aes(mutant, partner)) +
    geom_tile(aes(fill = deltarASA), color = 'white') +
    scale_fill_gradient(low = 'white', high = 'blue') +
    theme_grey() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = 'none') +
    ylab('delta rASA')
  
  # stack all plots onto one grid
  plot_list[[i]] <-
    plot_grid(emap_hm, emap_plot, GAP_plot, GEF_plot, GAPGEF_plot,
              apms_plot, iface_plot, nrow = 7, align = "v",
              rel_heights = c(2,3,1,1,1.5,2,2))
}

# set filename
if (GAP_GEF_complete_only & use_correlations) {
  fname <- 'images/corr_titrations_with_kinetics_and_apms_GAPGEF_Complete.pdf'
} else if (use_correlations) {
  fname <- 'images/corr_titrations_with_kinetics_and_apms.pdf'
} else if (GAP_GEF_complete_only) {
  fname <- 'images/score_titrations_with_kinetics_and_apms_GAPGEF_Complete.pdf'
} else {
  fname <- 'images/score_titrations_with_kinetics_and_apms.pdf'
}

pdf(fname, width = 10, height = 16)
for (p in plot_list) {print(p)}
dev.off()

