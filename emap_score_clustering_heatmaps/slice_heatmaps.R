#!/usr/bin/env Rscript --vanilla

clustfn_to_use <- commandArgs(trailingOnly=TRUE)[1]

library(tidyverse)
library(gplots)
library(gridGraphics)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
set.seed(94158)

make_Heatmap <- function(mat, title, name = 'legend',
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         row_dend_reorder = FALSE, column_dend_reorder = FALSE) {
  hm <- Heatmap(mat = mat,
                name=name,
                col = colorRamp2(c(-3, 0, 3), c("blue", "black", "yellow")),
                column_title = title,
                cluster_rows = cluster_rows,
                cluster_columns = cluster_columns,
                row_dend_reorder = row_dend_reorder,
                column_dend_reorder = column_dend_reorder)
  return(grid.grabExpr(draw(hm)))
}

filter_out_na <- function(.mat) {
  na_count <- apply(.mat, 1, function(x) sum(is.na(x)))
  .mat[which(ncol(.mat) - na_count >= 3),]
}


if (clustfn_to_use == 'pearson') {
  clustfn <- function(mat) {
    cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
    dissim <- as.dist((1 - cormat)/2)
    hc <- hclust(dissim, method = "average" )
    return(hc)
  }
}

if (clustfn_to_use == 'euclidean') {
  clustfn <- function(mat) {
    dissim <- dist(t(mat))
    hc <- hclust(dissim, method = "average" )
    return(hc)
  }
}
  
# get list of mutants for which we have GAP and GEF data
GG_mutants <-
  read_delim('datasets/GAP_GEF_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  arrange(GAP_GEF_ratio) %>% 
  pull(mutant)

# get ordering of mutants based on mass-spec
ms_order_GG <-
  scan('ms_ordering.txt', what=character(), quiet = TRUE) %>% 
  `[`(which(. %in% GG_mutants))


# save_heatmaps_to_pdf <- function(clustfn, pdf_name) {
# load and cluster emap, and subset of mutants with GAP/GEF data
emap <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols())

emap.mat <-
  emap %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
emap_GG.mat <- emap.mat[GG_mutants,]
emap.hc <- clustfn(t(emap.mat))
emap_GG.hc <- clustfn(t(emap_GG.mat))

# cluster interface ΔrASA
iface.mat <-
  read_delim('datasets/core_deltarASA_by_mutant_for_clustering.txt', delim = '\t', col_types = cols()) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()
iface_GG.mat <- iface.mat[GG_mutants[GG_mutants != 'GSP1-NAT'],]
iface.hc <- clustfn(t(iface.mat))
iface_GG.hc <- clustfn(t(iface_GG.mat))

# load gene sets
gene_sets <- read_delim('gene_sets_combined_min3.txt', delim = '\t', col_types = cols())
gene_set_list <- gene_sets %>% pull(gene_set) %>% unique()
plots_list <- list()


for (i in seq(length(gene_set_list))) {
    
  genes <-
    gene_sets %>%
    filter(gene_set == gene_set_list[i]) %>%
    filter(Gene %in% colnames(emap)) %>% 
    pull(Gene)
  
  print(paste(i, gene_set_list[i], length(genes), sep = '--')) 
  
  if (length(genes) < 3) {
    next
  }
  
  # load matrices with and without NAs, for all mutants and only those with GAP/GEF data
  mat.NA <- emap %>%
    select(mutant, genes) %>% 
    column_to_rownames(var="mutant") %>%
    as.matrix()
  mat_GG.NA <- mat.NA[GG_mutants,]
  
  mat <- filter_out_na(mat.NA)
  mat_GG <- filter_out_na(mat_GG.NA)
  
  # cluster
  row.hc <- clustfn(t(mat))
  col.hc <- clustfn(mat)
  row_GG.hc <- clustfn(t(mat_GG))
  col_GG.hc <- clustfn(mat_GG)
  
  # make a list to store heatmaps so they can be arranged on a grid
  
  gl <- list()
  
  # plot Heatmaps
  gl[[1]] <- make_Heatmap(mat = subset(emap.mat, select=colnames(mat.NA)), title = 'Full emap', name = 'S-score',
                          cluster_rows = emap.hc, cluster_columns = col.hc)
  gl[[2]] <- make_Heatmap(mat = mat, title = 'Recluster using only gene set', name = 'S-score',
                          cluster_rows = row.hc, cluster_columns = col.hc)
  gl[[3]] <- make_Heatmap(mat = mat.NA[!rownames(mat.NA) %in% c('GSP1-NAT','NTER3XFLAG WT','CTER3XFLAG WT'),],
                          title = 'Mutants ordered by interface', name = 'S-score',
                          cluster_rows = iface.hc, cluster_columns = col.hc)
  gl[[4]] <- make_Heatmap(mat = mat.NA, title = 'Seq order',name = 'S-score',
                          cluster_columns = col.hc)
  gl[[5]] <- grid.grabExpr(grid.text(paste0(gene_set_list[i], '\npearson correlation'), gp=gpar(fontsize=26)))
  gl[[6]] <- make_Heatmap(mat = subset(emap_GG.mat, select=colnames(mat_GG.NA)), title = 'Full emap', name = 'S-score',
                          cluster_rows = emap_GG.hc, cluster_columns = col.hc)
  gl[[7]] <- make_Heatmap(mat = mat_GG, title = 'Recluster using only gene set', name = 'S-score',
                          cluster_rows = row_GG.hc, cluster_columns = col_GG.hc)
  gl[[8]] <- make_Heatmap(mat = mat_GG.NA[!rownames(mat_GG.NA) %in% c('GSP1-NAT'),],
                          title = 'Mutants ordered by interface', name = 'S-score',
                          cluster_rows = iface_GG.hc, cluster_columns = col_GG.hc)
  gl[[9]] <- make_Heatmap(mat = mat_GG.NA, title = 'Mutants ordered by GAP/GEF efficiency ratio', name = 'S-score',
                          cluster_columns = col_GG.hc)
  gl[[10]] <- make_Heatmap(mat = mat_GG.NA[ms_order_GG,],
                          title = 'Mutants ordered by MS difference (GAP - GEF log2FC)', name = 'S-score',
                          cluster_rows = FALSE, cluster_columns = col_GG.hc)
  
  plots_list[[length(plots_list)+1]] <- arrangeGrob(grobs=gl, ncol=5)
}

pdfname <- paste0('images/heatmaps_', clustfn_to_use, '.pdf')

pdf(pdfname, height = 20, width = 30, onefile = TRUE)
for (p in plots_list) {
  grid.newpage()
  grid.draw(p)
}
graphics.off()
