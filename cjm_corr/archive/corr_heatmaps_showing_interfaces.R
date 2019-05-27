# version before May 23

##### load libraries
library(tidyverse)
library(gplots)
library(circlize)
library(dendextend)
library(psych)

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

##### load datafiles
load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/spitzemapko_correlations_and_pvalues_all.RData')
load('~/Box Sync/kortemmelab/home/tina/Gsp1/shared_datafiles/correlations/filtered_correlations.RData')
name2ORF <- read_delim('cjm_corr/spitzemap_name2ORF_index.txt', delim='\t')


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
  mutate('core' = case_when(interface == 'core' ~ T, interface != 'core' ~ F)) %>% 
  mutate('strain' = case_when(is.na(strain) ~ partner, T ~ strain)) %>% 
  select(strain, mutant, core) %>% 
  spread(strain, core) %>% 
  mutate('mog1' = NA) %>% 
  column_to_rownames('mutant') %>% 
  select(partner_index$strain) %>% 
  as.matrix()
core_res_mask[is.na(core_res_mask)] <- F

# make a matrix mask based on whether a point mutant is at a position
# whose ALA scan score was > 0.5 for a partner
alascan_mask <-
  read_tsv('cjm_corr/5-19-15tina_Gsp1_alascan_indexed.txt', col_types=cols()) %>% 
  select(yeast_index, yeast_wt_aa, protein, score, structure) %>% 
  left_join(partner_index, by = c('protein' = 'name')) %>%
  left_join(mutant_index, by = c('yeast_index' = 'yeastresnum')) %>% 
  filter(mutant %in% mutant_index$mutant) %>% 
  select(mutant, strain, score, structure) %>% 
  group_by(strain, mutant) %>% # YRB1 and KAP95 have two structures, so need to take max score
  mutate(max_score = max(score)) %>% # no negative scores less than -0.5
  ungroup() %>% 
  select(mutant, strain, max_score) %>% 
  unique() %>%
  spread(strain, max_score) %>% 
  mutate('mog1' = NA) %>% 
  column_to_rownames('mutant') %>% 
  select(-`<NA>`) %>% 
  as.matrix()
alascan_mask[is.na(alascan_mask)] <- 0
alascan_mask <- alascan_mask > 0.5

# Make a matrix of correlation values between partners and mutants
partners_vs_mutants <-
  correlations %>%
  filter(grepl('GSP1', query_uniq1)) %>% 
  filter(query_uniq2 %in% partner_index$strain) %>% 
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>% 
  filter(query_uniq1 %in% mutant_index$mutant)

mat <-
  partners_vs_mutants %>%
  select(query_uniq1, query_uniq2, pearson) %>% 
  spread(query_uniq2, pearson) %>% 
  column_to_rownames('query_uniq1') %>% 
  as.matrix()

# Make a pvalue mask based on a 0.05 threshold
# used to gray out heatmap cells for non-significant correlations
pmask <-
  partners_vs_mutants %>% 
  select(query_uniq1, query_uniq2, two_sided_p_value) %>% 
  spread(query_uniq2, two_sided_p_value) %>% 
  column_to_rownames('query_uniq1') %>% 
  as.matrix()
pmask <- pmask < 0.05 # so this masks for significant correlations
pmask[!pmask] <- NA
pmask <- 1*pmask

# cluster matrix of correlations and get orderings of rows (mutants) and columns (partners)
# cormat <- corr.test(t(mat), use = 'pairwise', method = 'pearson')
# dissim <- as.dist((1 - cormat$r)/2)
dissim <- dist(pmask*mat)
row_order <-
  cmdscale(dissim, eig = T, k = 1) %>% 
  `$`(point) %>% 
  as_tibble(rownames = 'mutant') %>% 
  arrange(V1) %>% 
  pull(mutant)
hc.row <- hclust(dissim, method = "single" )
hc.row <- rotate(hc.row, row_order)
# row_order <- hc.row$labels[hc.row$order]

# cols (partners)
# cormat <- corr.test(mat, use = 'pairwise', method = 'pearson')
# dissim <- as.dist((1 - cormat$r)/2)
dissim <- dist(t(pmask*mat))

col_order <-
  cmdscale(dissim, eig = T, k = 1) %>% 
  `$`(point) %>% 
  as_tibble(rownames = 'partner') %>% 
  arrange(V1) %>% 
  pull(partner)

hc.col <- hclust(dissim, method = "single")
hc.col <- rotate(hc.col, col_order)
# col_order <- hc.col$labels[hc.col$order]

# create a pdf to store all of the plots
# pdf('cjm_corr/mutant_partner_correlations.pdf', height = 10, width = 10)
# dev.off()
# all mutants, no filtering, boxes are core residues
title = 'Mutant-partner pearson corr,\nboxes = core residues (d_rASA)'
heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          # dendrogram = 'none', Rowv = F, Colv = F,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col)
          # add.expr={makeRects(core_res_mask[rev(row_order),col_order])}
          )

# all mutants, filtered, boxes are core residues
title = 'Mutant-partner pearson corr (pval<0.05),\nboxes = core residues (d_rASA)'
heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(core_res_mask[rev(row_order),col_order])})

# all mutants, no filtering, boxes are alanine scanning hits
title = 'Mutant-partner pearson corr,\nboxes =  score(ALA) > 0.5'
heatmap.2(mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(alascan_mask[rev(row_order),col_order])})

# all mutants, filtered, boxes are core residues
title = 'Mutant-partner pearson corr (pval<0.05),\nboxes =  score(ALA) > 0.5'
heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(alascan_mask[rev(row_order),col_order])})

# Remake plots, but only with strong mutants

# Make a matrix of correlation values between partners and mutants

partners_vs_mutants <-
  correlations %>%
  filter(grepl('GSP1', query_uniq1)) %>% 
  filter(query_uniq2 %in% partner_index$strain) %>% 
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>% 
  filter(query_uniq1 %in% filter(mutant_index, is_strong)$mutant)

mat <-
  partners_vs_mutants %>%
  select(query_uniq1, query_uniq2, pearson) %>% 
  spread(query_uniq2, pearson) %>% 
  column_to_rownames('query_uniq1') %>% 
  as.matrix()

# Make a pvalue mask based on a 0.05 threshold
# used to gray out heatmap cells for non-significant correlations
pmask <-
  partners_vs_mutants %>% 
  select(query_uniq1, query_uniq2, two_sided_p_value) %>% 
  spread(query_uniq2, two_sided_p_value) %>% 
  column_to_rownames('query_uniq1') %>% 
  as.matrix()
pmask <- pmask < 0.05 # so this masks for significant correlations
pmask[!pmask] <- NA
pmask <- 1*pmask

# cluster matrix of correlations and get orderings of rows (mutants) and columns (partners)
cormat <- corr.test(t(mat), use = 'pairwise', method = 'pearson')
dissim <- as.dist((1 - cormat$r)/2)
hc.row <- hclust(dissim, method = "single" )
row_order <- hc.row$labels[hc.row$order]

# cols (partners)
cormat <- corr.test(mat, use = 'pairwise', method = 'pearson')
dissim <- as.dist((1 - cormat$r)/2)
hc.col <- hclust(dissim, method = "single" )
col_order <- hc.col$labels[hc.col$order]

# all mutants, no filtering, boxes are core residues
title = 'Mutant-partner pearson corr,\nstrong mutants only\nboxes = core residues (d_rASA)'
heatmap.2(mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(core_res_mask[rev(row_order),col_order])})

# all mutants, filtered, boxes are core residues
title = 'Mutant-partner pearson corr (pval<0.05),\nstrong mutants only\nboxes = core residues (d_rASA)'
heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(core_res_mask[rev(row_order),col_order])})

# all mutants, no filtering, boxes are alanine scanning hits
title = 'Mutant-partner pearson corr,\nstrong mutants only\nboxes =  score(ALA) > 0.5'
heatmap.2(mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(alascan_mask[rev(row_order),col_order])})

# all mutants, filtered, boxes are core residues
title = 'Mutant-partner pearson corr (pval<0.05),\nstrong mutants only\nboxes =  score(ALA) > 0.5'
heatmap.2(pmask*mat, trace = 'none', col = colramp, na.color='gray', srtCol=45,
          lhei=c(1,4), lwid=c(1,2), main = title,
          Rowv = as.dendrogram(hc.row), Colv = as.dendrogram(hc.col),
          add.expr={makeRects(alascan_mask[rev(row_order),col_order])})
dev.off()
