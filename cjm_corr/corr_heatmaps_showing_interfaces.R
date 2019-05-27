
##### load libraries
library(tidyverse)
library(gplots)
library(circlize)
library(dendextend)
library(psych)
library(ggforce)
library(RColorBrewer)

library(ggpubr)




# color for correlations
colramp = colorRampPalette(c("red", "white", "green"))

# Function for making selection rectangles around selection cells
makeRects <- function(mask){
  coords = expand.grid(dim(mask)[1]:1, 1:dim(mask)[2])[mask,]
  xl=coords[,2]-0.49
  yb=coords[,1]-0.49
  xr=coords[,2]+0.49
  yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=3)
}

# Function to make heatmap with boxes
make_heatmaps <- function(table) {

  mat <-
    table %>%
    select(query_uniq1, query_uniq2, pearson) %>% 
    spread(query_uniq2, pearson) %>% 
    column_to_rownames('query_uniq1') %>% 
    as.matrix()
  
  # Make a pvalue mask based on a 0.05 threshold
  # used to gray out heatmap cells for non-significant correlations
  pmask <-
    table %>% 
    select(query_uniq1, query_uniq2, two_sided_p_value) %>% 
    spread(query_uniq2, two_sided_p_value) %>% 
    column_to_rownames('query_uniq1') %>% 
    as.matrix()
  pmask <- pmask < 0.05 # so this masks for significant correlations
  pmask[!pmask] <- NA
  pmask <- 1*pmask
  
  # cluster matrix of correlations and get orderings of rows (mutants) and columns (partners)
  cor.row <- corr.test(t(mat), use = 'pairwise', method = 'pearson')
  dissim <- as.dist((1 - cor.row$r)/2)
  hc.row <- hclust(dissim, method = "single" )
  
  cor.col <- corr.test(mat, use = 'pairwise', method = 'pearson')
  dissim <- as.dist((1 - cor.col$r)/2)
  hc.col <- hclust(dissim, method = "single" )
  
  # use the first principal coordinate (classical multidimensional scaling) to flip leaves
  row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'mutant') %>%
    arrange(`V1`) %>% pull(mutant)
  hc.row <- rotate(hc.row, row_pcoor)
  
  col_pcoor <-
    cmdscale(dist(t(mat)), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'mutant') %>%
    arrange(`V1`) %>% pull(mutant)
  hc.col <- rotate(hc.col, col_pcoor)
  
  # get orderings (needed to place boxes)
  row_order <- hc.row$labels[hc.row$order]
  col_order <- hc.col$labels[hc.col$order]
  
  # all mutants, no filtering, boxes are core residues
  title = 'Mutant-partner pearson corr,\nboxes = core residues (d_rASA)'
  heatmap.2(mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
            lhei=c(1,4), lwid=c(1,2), main = title,
            Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
            add.expr={makeRects(core_res_mask[rev(row_order),col_order])})
  
  # all mutants, filtered, boxes are core residues
  title = 'Mutant-partner pearson corr (pval<0.05),\nboxes = core residues (d_rASA)'
  heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
            lhei=c(1,4), lwid=c(1,2), main = title,
            Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
            add.expr={makeRects(core_res_mask[rev(row_order),col_order])})
}


##### load datafiles
load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/spitzemapko_correlations_and_pvalues_all.RData')
load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')
name2ORF <- read_delim('spitzemap_name2ORF_index.txt', delim='\t') # in cjm_corr

##### Make index files for Gsp1 mutants and partners
mutant_index <-
  filter(name2ORF, grepl('GSP1', name)) %>% 
  mutate('is_strong' = case_when(name %in% filtered_correlations$query_uniq1 ~ T, T ~ F)) %>% 
  rename('strain' = name) %>% 
  separate(strain, sep = ' - ', into = c('GSP1','mutant')) %>% 
  select(-GSP1, -ORF) %>% 
  filter(! mutant %in% c('GSP1-NAT','NTER3XFLAG WT', 'CTER3XFLAG WT')) %>% 
  mutate('yeastresnum' = as.numeric(str_sub(mutant, 2, nchar(mutant)-1)))
  
partner_names <- c('MSN5','SRP1','LOS1','YRB1','YRB2','KAP95','RNA1',
                   'SRM1','MTR10','PSE1','NTF2','CRM1','CSE1','KAP104','MOG1')
partner_ORFs <- c('YDR335W','YNL189W','YKL205W','YDR002W','YIL063C',
                  'YLR347C','YMR235C','YGL097W','YOR160W','YMR308C',
                  'YER009W','YGR218W','YGL238W','YBR017C','YJR074W') 
partner_index <- data.frame(partner_names, partner_ORFs)
names(partner_index) <- c('name', 'ORF')
partner_strains <- filter(name2ORF, str_detect(ORF, paste(partner_ORFs, collapse = '|')))
partner_index <-
  left_join(partner_index, partner_strains, by = 'ORF') %>% 
  rename('name' = 'name.x', 'strain' = 'name.y') %>% 
  filter(strain %in% correlations$query_uniq1)

# make a matrix mask based on whether a point mutant is at an interface
# core position for a partner
core_res_mask <-
  read_tsv('~/Box Sync/kortemmelab/home/cjmathy_lab/gsp1/ran_structures/SASA_interfaces.txt', col_types = cols()) %>% 
  left_join(partner_index, by = c('partner' = 'name')) %>% 
  left_join(mutant_index, by = 'yeastresnum') %>% 
  filter(mutant %in% mutant_index$mutant) %>% 
  mutate('is_core' = case_when(interface == 'core' ~ T, interface != 'core' ~ F)) %>% 
  mutate('strain' = case_when(is.na(strain) ~ partner, T ~ strain)) %>% 
  select(strain, mutant, is_core) %>% 
  spread(strain, is_core) %>% 
  mutate('mog1' = NA) %>% 
  column_to_rownames('mutant') %>% 
  select(partner_index$strain) %>% 
  as.matrix()
core_res_mask[is.na(core_res_mask)] <- F

# Make a table of correlation values between partners and mutants
partners_vs_all_mutants <-
  correlations %>%
  filter(grepl('GSP1', query_uniq1)) %>% 
  filter(query_uniq2 %in% partner_index$strain) %>% 
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>% 
  filter(query_uniq1 %in% mutant_index$mutant)

# Also make a table for just strong mutants
partners_vs_strong_mutants <-
  correlations %>%
  filter(grepl('GSP1', query_uniq1)) %>% 
  filter(query_uniq2 %in% partner_index$strain) %>% 
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>% 
  filter(query_uniq1 %in% filter(mutant_index, is_strong)$mutant)

# create a pdf to store all of the plots
pdf('cjm_corr/mutant_partner_correlations.pdf', height = 10, width = 10)

# make heatmaps with and without pval masking, for all mutants and strong mutants
make_heatmaps(partners_vs_all_mutants)
make_heatmaps(partners_vs_strong_mutants)

dev.off()

# Sina plots, in interface, not in interface
core_res_table <-
  read_tsv('~/Box Sync/kortemmelab/home/cjmathy_lab/gsp1/ran_structures/SASA_interfaces.txt', col_types = cols()) %>%
  left_join(partner_index, by = c('partner' = 'name')) %>% 
  left_join(mutant_index, by = 'yeastresnum') %>% 
  filter(mutant %in% mutant_index$mutant) %>% 
  mutate('is_core' = case_when(interface == 'core' ~ T, interface != 'core' ~ F)) %>% 
  mutate('strain' = case_when(is.na(strain) ~ partner, T ~ strain)) %>% 
  filter(is_core) %>% 
  select(strain, mutant, is_core)
  
strong_mutants <-
  mutant_index %>% 
  select(mutant, is_strong)

pdf('cjm_corr/sinaplots.pdf', height = 6, width = 8)

sinaplot_table <-
  partners_vs_all_mutants %>%
  select(query_uniq1, query_uniq2, pearson) %>% 
  rename('mutant' = query_uniq1, 'strain' = query_uniq2) %>% 
  left_join(core_res_table) %>%
  left_join(strong_mutants) %>% 
  mutate(is_core = case_when(is_core ~ 'Mutant in interface',
                             is.na(is_core) ~ 'Mutant not in interface')) %>% 
  mutate(is_core = factor(is_core, levels = c('Mutant not in interface', 'Mutant in interface'))) %>% 
  rename('partner' = strain)

title1 = paste0('Correlations between all mutants and core partners',
                '\nResidues deemed in interface based on d_rASA')

title2 = paste0('Correlations between strong mutants and core partners',
                '\nResidues deemed in interface based on d_rASA')


colors1 = brewer.pal(n = length(unique(sinaplot_table$partner)), name = 'Set1')
colors2 = brewer.pal(n = length(unique(sinaplot_table$partner)), name = 'Set2')
palette = c(colors1, colors2)

pdf('sina_plot.pdf')
set.seed(3)
ggplot(filter(sinaplot_table, is_strong), aes(x = is_core, y = pearson, color = is_core)) +
  geom_sina(alpha = 0.8, size = 3, maxwidth = 0.8) +
  scale_color_manual(name = 'group', values = c('black', 'salmon'),
                   labels = c('weak mutants','strong mutants')) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), color = 'red',
               geom="errorbar", width=0.2, size = 2) +
  stat_summary(fun.y=mean, geom="point", color="red", size = 4) +
  ylim(c(-0.15,0.5)) + xlab('') + ylab('Pearson correlation, r') + ggtitle(title1) +
  theme(text=element_text(size=40,  family="Helvetica"),
        ) %>% print()
dev.off()





pdf('dotplots.pdf')

# all mutants
ggplot(sinaplot_table, aes(x = is_core, y = pearson)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.02) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), color = 'red',
               geom="errorbar", width=0.2, size = 2) +
  stat_summary(fun.y=mean, geom="point", color="red", size = 4) +
  stat_compare_means(label.x = 1.3, label.y = 0.5, size = 6) +
  ylim(c(-0.3,0.5)) + xlab('') + ylab('Pearson correlation, r') + ggtitle(title1) +
  theme(text=element_text(size=16,  family="Helvetica")) %>% print()


# Just strong mutants
ggplot(filter(sinaplot_table, is_strong), aes(x = is_core, y = pearson)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.02) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), color = 'red',
               geom="errorbar", width=0.2, size = 2) +
  stat_summary(fun.y=mean, geom="point", color="red", size = 4) +
  stat_compare_means(label.x = 1.3, label.y = 0.5, size = 6) +
  ylim(c(-0.3,0.5)) + xlab('') + ylab('Pearson correlation, r') + ggtitle(title2) +
  theme(text=element_text(size=16,  family="Helvetica")) %>% print()

dev.off()

ggplot(sinaplot_table) +
  geom_dotplot(aes(x = is_core, y = pearson, fill = is_strong),
               binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.02) +
  stat_summary(aes(x = is_core, y = pearson), fun.data=mean_sdl, fun.args = list(mult=1), color = 'red',
               geom="errorbar", width=0.2, size = 2) +
  # stat_summary(fun.y=mean, geom="point", color="red", size = 4) +
  xlab('') + ylab('Pearson correlation, r') + ggtitle(title) +
  scale_fill_manual(name = 'group', values = palette,
                    labels = c('weak mutants','strong mutants')) +
  theme(text=element_text(size=16,  family="Helvetica")) %>% print()

ggplot(sinaplot_table, aes(x = is_core, y = pearson, color = partner)) +
  geom_sina(size = 2) + xlab('') + ylab('Pearson correlation, r') + ggtitle(title) +
  scale_color_manual(values = palette) +
  theme(text=element_text(size=16,  family="Helvetica")) %>% print()

ggplot(sinaplot_table, aes(x = is_core, y = pearson, shape = is_strong, color = partner)) +
  geom_sina(size = 2) + xlab('') + ylab('Pearson correlation, r') + ggtitle(title) +
  scale_color_manual(values = palette) +
  scale_shape_manual(name = 'group',
                     labels = c('weak mutants','strong mutants'),
                     values = c('circle','triangle')) +
  theme(text=element_text(size=16,  family="Helvetica")) %>% print()

dev.off()



# correlations of correlations, mutants and partners


colramp = colorRampPalette(c("orange", "white", "purple"))




GSP1_names <- filter(name2ORF, grepl('GSP1', name)) %>% pull(name)
core_partners <- c('rna1-s116f','rna1-1','srp1-5001',
                   'yrb1-51','kap95-e126k','srm1-g282s',
                   'srm1-ts','ntf2-5001','ntf2-h104y','mog1')




# compute correlations of correlations matrix, filtering correlations by pval < 0.05 first
corr_of_corr <-
  filtered_correlations %>%
  filter(two_sided_p_value < 0.05) %>%
  select(query_uniq1, query_uniq2, pearson) %>%
  spread(query_uniq2, pearson) %>%
  column_to_rownames('query_uniq1') %>%
  as.matrix() %>%
  t() %>%
  corr.test(., use = "pairwise", method = "pearson", ci=F)

# make mask for corr of corr's with pval < 0.05
pmask <- corr_of_corr$p 
pmask[upper.tri(pmask)] = 0 # use the unadjusted p-values
pmask = pmask + t(pmask)
pmask <- 1*(pmask < 0.05) # threshold for p< 0.05
pmask[pmask==0] <- NA

# get smaller matrices, of just Gsp1 mutants and partners
strong_mutants <-
  mutant_index %>% 
  select(mutant, is_strong)

mutants <- intersect(rownames(corr_of_corr$r), GSP1_names)
mutants <- setdiff(mutants, c('GSP1 - CTER3XFLAG WT','GSP1 - NTER3XFLAG WT'))


mat <- corr_of_corr$r[mutants, core_partners]
rownames(mat) <- lapply(rownames(mat), function(x) substr(x, 8, nchar(x)))
pmask <- pmask[mutants, core_partners]
rownames(pmask) <- lapply(rownames(pmask), function(x) substr(x, 8, nchar(x)))

# cluster matrix of corr of corrs using euclidean (small number of partners) and get orderings of rows (mutants) and columns (partners)

# rows (mutants)
dissim <- dist(pmask*mat)
hc.row <- hclust(dissim, method = "single" )
row_order <- hc.row$labels[hc.row$order]

# cols (partners)
dissim <- dist(t(pmask*mat))
hc.col <- hclust(dissim, method = "single" )
col_order <- hc.col$labels[hc.col$order]

heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = 'title',
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          # Colv = core_partners,
          # dendrogram='row'
          # add.expr={makeRects(core_res_mask[rev(row_order),col_order])}
)

dev.off()



pdf('cjm_corr/corr_of_corr_partners_mutants.pdf', height = 7, width = 7)

Heatmap(cormat_filtered[GSP1_to_keep, core_partners],
        name = 'r',
        col = colramp,
        column_title = 'Correlations of correlations',
        cluster_columns = FALSE,
        column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
        row_dend_width = unit(30, "mm"), column_dend_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

dev.off()
