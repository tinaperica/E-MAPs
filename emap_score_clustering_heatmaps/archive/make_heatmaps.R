library(tidyverse)  # data manipulation
library(gplots)
library(gridGraphics)
library(grid)
library(gridExtra)
library(RColorBrewer)
# library(cluster)    # clustering algorithms
# library(factoextra) # clustering visualization
# library(dendextend) # for comparing two dendrograms
set.seed(94158)
setwd('/Users/cjmathy/lab/gsp1/emaps/raw_gsp1_emap_clustering')

logfile <- file('logfile.txt', open='wt')
sink(logfile, type='message')



# prepare list of mutant orderings for displaying heatmaps in different row orders
mut_orderings <- list()
mut_orderings$seq_order <- scan('gsp1_mutants_ordered_by_residue.txt', what="", sep="\n")

# get list of mutants for which we have GAP and GEF data
mutants_with_GAP_GEF_data <- read_delim('datasets/GAP_GEF_ratio_for_clustering.txt', delim = '\t', col_types = cols())$mutant

# define clustering function
clustfn <- function(mat) {
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  hc <- hclust(dissim, method = "average" )
  return(hc)
}


#--cluster full emap

# Read in full emap for clustering
dset <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols())

# cluster mutants
mut.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
colnames(mut.mat) <- dset %>% pull(mutant)
mut.hc <- clustfn(mut.mat)

# save mutant orderings
ordering <- mut.hc$labels[mut.hc$order]
mut_orderings$full_emap <- ordering
mut_orderings$full_emap_GAPGEF_cluster_then_filter <- ordering[ordering %in% mutants_with_GAP_GEF_data]

# cluster arrays
arr.mat <- select(dset, -mutant) %>% as.matrix()
rownames(arr.mat) <- dset %>% pull(mutant)
arr.hc <- clustfn(arr.mat)

# make a heatmap of the full emap, with mutants and arrays clustered
# pdf('images/gsp1_emap_full.pdf', height = 10, width = 60)
# heatmap.2(arr.mat,
#           Rowv = as.dendrogram(mut.hc), Colv = as.dendrogram(arr.hc),
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101),
#           trace='none', density.info="none",
#           lhei = c(1.5, 6), lwid = c(0.25, 6),
#           cexCol = 0.3, cexRow = 1)
# dev.off()

# also plot emap without mutant clustering (previously ordered by residue)
# pdf('images/gsp1_emap_full_mutants_ordered_by_res.pdf', height = 10, width = 60)
# heatmap.2(arr.mat,
#           Rowv = FALSE, Colv = as.dendrogram(arr.hc),
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101),
#           trace='none', density.info="none",
#           lhei = c(1.5, 6), lwid = c(0.25, 6),
#           cexCol = 0.3, cexRow = 1)
# dev.off()

# store trees in better described objects
full_emap_rowtree <- mut.hc
full_emap_coltree <- arr.hc

#--cluster full emap, GAPGEF only
# Read in full emap for clustering
dset <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
  filter(mutant %in% mutants_with_GAP_GEF_data)

# cluster mutants
mut.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
colnames(mut.mat) <- dset %>% pull(mutant)
mut.hc <- clustfn(mut.mat)

# save mutant orderings
ordering <- mut.hc$labels[mut.hc$order]
mut_orderings$full_emap_GAPGEF_filter_then_cluster <- ordering

# cluster arrays
arr.mat <- select(dset, -mutant) %>% as.matrix()
rownames(arr.mat) <- dset %>% pull(mutant)

# filtering for GAP/GEF completeness leaves 3 pairs of array gene
# correlations as NA, so replace these with maximum distance values (1.0)
cormat <- cor(arr.mat, use = "pairwise.complete.obs", method = "pearson")
dissim <- as.dist((1 - cormat)/2)
dissim[which(sapply(dissim, function(x)all(is.na(x))))] <- 1.0
arr.hc <- hclust(dissim, method = "average" )

# make a heatmap of the full emap, with mutants and arrays clustered
# pdf('images/gsp1_emap_full_GAPGEF_only.pdf', height = 10, width = 60)
# heatmap.2(arr.mat,
#           Rowv = as.dendrogram(mut.hc), Colv = as.dendrogram(arr.hc),
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101),
#           trace='none', density.info="none",
#           lhei = c(1.5, 6), lwid = c(0.25, 6),
#           cexCol = 0.3, cexRow = 1)
# dev.off()

# # also plot emap without mutant clustering (previously ordered by residue)
# pdf('images/gsp1_emap_full_GAPGEF_only_mutants_ordered_by_res.pdf', height = 10, width = 60)
# heatmap.2(arr.mat,
#           Rowv = FALSE, Colv = as.dendrogram(arr.hc),
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101),
#           trace='none', density.info="none",
#           lhei = c(1.5, 6), lwid = c(0.25, 6),
#           cexCol = 0.3, cexRow = 1)
# dev.off()

# store trees in better described objects
full_emap_GAP_GEF_rowtree <- mut.hc
full_emap_GAP_GEF_coltree <- arr.hc


#--cluster interface dataset


# Read in full emap for clustering
dset <- read_delim('datasets/core_deltarASA_by_mutant_for_clustering.txt', delim = '\t', col_types = cols())

# cluster mutants
mut.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
colnames(mut.mat) <- dset %>% pull(mutant)
cormat <- cor(mut.mat, use = "everything", method = "pearson")
cormat[1:3,] = 1
cormat[,1:3] = 1
dissim <- as.dist((1 - cormat)/2)
mut.hc <- hclust(dissim, method = "average" )

# save mutant orderings
ordering <- mut.hc$labels[mut.hc$order]
mut_orderings$interface <- ordering
mut_orderings$interface_GAPGEF_cluster_then_filter <- ordering[ordering %in% mutants_with_GAP_GEF_data]

# cluster arrays
partner <- list()
partner.mat <- select(dset, -mutant) %>% as.matrix()
rownames(partner.mat) <- dset %>% pull(mutant)
partner.hc <- clustfn(partner.mat)

# heatmap
# pdf('images/interface_clustering.pdf', height = 10, width = 8)
# heatmap.2(partner.mat,
#           Rowv = as.dendrogram(mut.hc), Colv = as.dendrogram(partner.hc),
#           breaks = seq(0, 1, length.out = 101),
#           col = colorRampPalette(c("black","cyan")),
#           trace='none', density.info="none",
#           main = 'Clustered by change in ASA after binding partner',
#           lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0)),
#           lhei = c(0.5,4,0.75), lwid = c(0.5,4,0.25),
#           cexCol = 1.5, cexRow = 1.25, srtCol=45)
# dev.off()

# store trees in better described objects
interface_rowtree <- mut.hc
interface_coltree <- partner.hc

#--cluster interface dataset, GAPGEF only
# Read in full emap for clustering
dset <- read_delim('datasets/core_deltarASA_by_mutant_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
  filter(mutant %in% mutants_with_GAP_GEF_data)

# cluster mutants
mut.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
colnames(mut.mat) <- dset %>% pull(mutant)
cormat <- cor(mut.mat, use = "everything", method = "pearson")
cormat[1:3,] = 1
cormat[,1:3] = 1
dissim <- as.dist((1 - cormat)/2)
mut.hc <- hclust(dissim, method = "average" )

# save mutant orderings
ordering <- mut.hc$labels[mut.hc$order]
mut_orderings$interface_GAPGEF_filter_then_cluster <- ordering

# cluster arrays
partner <- list()
partner.mat <- select(dset, -mutant) %>% as.matrix()
rownames(partner.mat) <- dset %>% pull(mutant)
partner.hc <- clustfn(partner.mat)

# heatmap
# pdf('images/interface_GAPGEF_only_clustering.pdf', height = 10, width = 8)
# heatmap.2(partner.mat,
#           Rowv = as.dendrogram(mut.hc), Colv = as.dendrogram(partner.hc),
#           breaks = seq(0, 1, length.out = 101),
#           col = colorRampPalette(c("black","cyan")),
#           trace='none', density.info="none",
#           main = 'Clustered by change in ASA after binding partner',
#           lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0)),
#           lhei = c(0.5,4,0.75), lwid = c(0.5,4,0.25),
#           cexCol = 1.5, cexRow = 1.25, srtCol=45)
# dev.off()

# store trees in better described objects
interface_GAP_GEF_rowtree <- mut.hc
interface_GAP_GEF_coltree <- partner.hc


#--cluster GAP/GEF dataset
# Read in full emap for clustering
dset <- read_delim('datasets/GAP_GEF_effic_for_clustering.txt', delim = '\t', col_types = cols())

# cluster mutants
mut.mat <- select(dset, -mutant) %>% as.matrix()
rownames(mut.mat) <- dset %>% pull(mutant)
distmat <- dist(mut.mat, method = "euclidean")
mut.hc <- hclust(distmat, method = "average")

mut_orderings$GAPGEF_effic <- mut.hc$labels[mut.hc$order]

# make a heatmap of the full emap, with mutants and arrays clustered
# pdf('images/GAPGEF_clustering.pdf', height = 10, width = 6)
# heatmap.2(mut.mat,
#           Rowv = as.dendrogram(mut.hc), Colv = FALSE,
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101),
#           trace='none', density.info="none",
#           lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0)),
#           lhei = c(0.5,4,1), lwid = c(0.5,2,0.25),
#           cexCol = 1, cexRow = 1, srtCol=0, adjCol = 0.5)
# dev.off()

GAP_GEF_rowtree <- mut.hc


# Cut EMAP into k clusters, and consider only those clusters with at least 3
# members. After examining the SGD descriptions, about half (21/44) strongly
# enrich one biological process or complex. These clusters were manually labeled
# and the file 'subset_cluster_names.txt' was made.
# NOTE: different clusters if you filter to GAP/GEF complete first

# set k = number of clusters to make
# k = 1200
# cluster_assignements <- cutree(full_emap_coltree, k = k)
# # cluster_assignements <- cutree(full_emap_GAP_GEF_coltree, k = k)
# clusters_with_descriptions <-
#   data.frame('name' = names(cluster_assignements),
#              'cluster_number' = unname(cluster_assignements),
#              stringsAsFactors = FALSE) %>% 
#   left_join(read_delim('SGD_descriptions.txt', delim='\t', col_types = cols()), by = 'name')

# clusters_with_descriptions %>% 
#   write_delim('unbiased_array_clusters_with_descriptions.txt', delim='\t')

# cluster_dir = 'clusters/'
# dir.create(file.path(getwd(), cluster_dir))

# clusters_with_descriptions %>% 
#   group_by(cluster_number) %>% 
#   filter(n() > 2) %>% 
#   nest(.key = 'cluster_data') %>% 
#   mutate(data = walk2(cluster_data,
#                       paste0(cluster_dir, '/', cluster_number, '.csv'),
#                       delim = ',', na = '', write_delim))


#-- Make gene sets
# gene_sets_from_clustering <- 
#   clusters_with_descriptions %>% 
#   inner_join(read_delim('subset_cluster_names.txt', delim = ',', col_types = cols())) %>% 
#   arrange(cluster_number) %>% 
#   select(ORF, name, cluster_name) %>% 
#   rename('Gene' = name, 'gene_set'  = cluster_name) %>% 
#   mutate('set_label' = 'from_clustering')

# gene_sets_processes <-
#   read_delim('gene_sets_for_gsp1_functions.txt', delim='\t', col_types = cols()) %>% 
#   rename('gene_set' = term) %>% 
#   mutate('set_label' = 'curated_processes')

# gene_sets_complexes <-
#   read_delim('all_wodak_complexes.txt', delim='\t', col_types = cols()) %>% 
#   rename('Gene' = Name, 'gene_set' = Complex) %>% 
#   mutate('set_label' = 'wodak_complexes')

# gene_sets <- bind_rows(gene_sets_from_clustering,
#                        gene_sets_processes,
#                        gene_sets_complexes)

# write_delim(gene_sets, 'gene_sets_combined.txt', delim='\t')
gene_sets <- read_delim('gene_sets_combined.txt', delim='\t', col_types = cols())


#--Plot heatmaps, euclidean
make_heatmap <- function(mat, Rowv, Colv, main) {
  heatmap.2(mat, Rowv = Rowv, Colv = Colv, main = main,
          col = colorRampPalette(c('blue','black','yellow')),
          breaks = seq(-3, 3, length.out = 101), trace='none', density.info="none", keysize = 0.2,
          lmat = rbind(c(4,3,0),c(2,1,0),c(0,0,0)), lhei = c(2, 5, 0.1), lwid = c(1, 2, 0.5),
          cexCol = 1, cexRow = 1, srtCol=45)
  grid.echo()
  return(grid.grab())
}

gene_set_list <- gene_sets %>% pull(gene_set) %>% unique()
plots_list <- list()

for (i in seq(length(gene_set_list))) {

  genes <-
    read_delim('gene_sets_combined.txt', delim = '\t', col_types = cols()) %>% 
    filter(gene_set == gene_set_list[i]) %>%
    pull(Gene)

  # complete matrices, for plotting comparison dendrograms
  dset.withNA <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
    gather(key = 'gene', value = score, -mutant) %>%
    filter(gene %in% genes) %>%
    spread(key = gene, value = score)
  mut.mat.withNA <- select(dset.withNA, -mutant) %>% as.matrix()
  rownames(mut.mat.withNA) <- dset.withNA %>% pull(mutant)
  dset_GG.withNA <- dset.withNA %>% filter(mutant %in% mutants_with_GAP_GEF_data)
  mut_GG.mat.withNA <- select(dset_GG.withNA, -mutant) %>% as.matrix()
  rownames(mut_GG.mat.withNA) <- dset_GG.withNA %>% pull(mutant)
  
  # make dataset without NAs for euclidean distance
  dset <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
    gather(key = 'gene', value = score, -mutant) %>%
    filter(gene %in% genes) %>%
    spread(key = gene, value = score) %>% 
    drop_na() # need to get rid of rows with all NA
  
  mut.mat <- select(dset, -mutant) %>% as.matrix()
  rownames(mut.mat) <- dset %>% pull(mutant)
  mut.dist <- dist(mut.mat, method = "euclidean")
  mut.hc <- hclust(mut.dist, method = "average")
  arr.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
  colnames(arr.mat) <- dset %>% pull(mutant)
  arr.dist <- dist(arr.mat, method = "euclidean")
  arr.hc <- hclust(arr.dist, method = "average")
 
  # GAP GEF dataset
  dset_GG <- dset %>% filter(mutant %in% mutants_with_GAP_GEF_data)
  mut_GG.mat <- select(dset_GG, -mutant) %>% as.matrix()
  rownames(mut_GG.mat) <- dset_GG %>% pull(mutant)
  mut_GG.dist <- dist(mut_GG.mat, method = "euclidean")
  mut_GG.hc <- hclust(mut_GG.dist, method = "average")
  arr_GG.mat <- select(dset_GG, -mutant) %>% as.matrix() %>% t()
  colnames(arr_GG.mat) <- dset_GG %>% pull(mutant)
  arr_GG.dist <- dist(arr_GG.mat, method = "euclidean")
  arr_GG.hc <- hclust(arr_GG.dist, method = "average")


  gl <- gList()
  gl[[1]] <- make_heatmap(mat = mut.mat.withNA, main = 'Full emap',
                          Rowv = as.dendrogram(full_emap_rowtree),
                          Colv = as.dendrogram(arr.hc))
  break
  warnings()

  gl[[2]] <- make_heatmap(mat = mut.mat.withNA, main = 'Seq order',
                          Rowv = FALSE,
                          Colv = as.dendrogram(arr.hc))
  warnings()
  gl[[3]] <- make_heatmap(mat = mut.mat, main = 'Gene set recluster',
                          Rowv = as.dendrogram(mut.hc),
                          Colv = as.dendrogram(arr.hc))
  warnings()
  gl[[4]] <- make_heatmap(mat = mut.mat.withNA, main = 'Interface',
                          Rowv = as.dendrogram(interface_rowtree),
                          Colv = as.dendrogram(arr.hc))
  warnings()
  # gl[[1]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
  # gl[[2]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
  # gl[[3]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
  # gl[[4]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
  gl[[5]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
  gl[[6]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Full emap',
                          Rowv = as.dendrogram(full_emap_GAP_GEF_rowtree),
                          Colv = as.dendrogram(arr_GG.hc))
  gl[[7]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Seq order',
                          Rowv = FALSE,
                          Colv = as.dendrogram(arr_GG.hc))
  gl[[8]] <- make_heatmap(mat = mut_GG.mat, main = 'Gene set recluster',
                          Rowv = as.dendrogram(mut_GG.hc),
                          Colv = as.dendrogram(arr_GG.hc))
  gl[[9]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Interface',
                          Rowv = as.dendrogram(interface_GAP_GEF_rowtree),
                          Colv = as.dendrogram(arr_GG.hc))
  gl[[10]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'GAP GEF',
                          Rowv = as.dendrogram(GAP_GEF_rowtree),
                          Colv = as.dendrogram(arr_GG.hc))
  if(i == 2){break}

  grid.arrange(grobs=gl, ncol=5, clip=TRUE, newpage=TRUE)
  plots_list[[i]] <- recordPlot()
}
pdf('images/heatmaps_euclid.pdf', height = 20, width = 30, onefile = TRUE)
for (my.plot in plots_list) {
  replayPlot(my.plot)
}
graphics.off()
warnings()
sink(type='message')
close(logfile)

# #--Plot heatmaps, correlation
# make_heatmap <- function(mat, Rowv, Colv, main) {
#   heatmap.2(mat, Rowv = Rowv, Colv = Colv, main = main,
#           col = colorRampPalette(c('blue','black','yellow')),
#           breaks = seq(-3, 3, length.out = 101), trace='none', density.info="none", keysize = 0.2,
#           lmat = rbind(c(4,3,0),c(2,1,0),c(0,0,0)), lhei = c(2, 5, 0.1), lwid = c(1, 2, 0.5),
#           cexCol = 1, cexRow = 1, srtCol=45)
#   grid.echo()
#   return(grid.grab())
# }

# gene_set_list <- gene_sets %>% pull(gene_set) %>% unique()
# plots_list <- list()
# for (i in seq(length(gene_set_list))) {

#   genes <-
#     read_delim('gene_sets_combined.txt', delim = '\t', col_types = cols()) %>% 
#     filter(gene_set == gene_set_list[i]) %>%
#     pull(Gene)

#   # complete matrices, for plotting comparison dendrograms
#   dset.withNA <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
#     gather(key = 'gene', value = score, -mutant) %>%
#     filter(gene %in% genes) %>%
#     spread(key = gene, value = score)
#   mut.mat.withNA <- select(dset.withNA, -mutant) %>% as.matrix()
#   rownames(mut.mat.withNA) <- dset.withNA %>% pull(mutant)
#   dset_GG.withNA <- dset.withNA %>% filter(mutant %in% mutants_with_GAP_GEF_data)
#   mut_GG.mat.withNA <- select(dset_GG.withNA, -mutant) %>% as.matrix()
#   rownames(mut_GG.mat.withNA) <- dset_GG.withNA %>% pull(mutant)
  
#   # make dataset without NAs for euclidean distance
#   dset <- read_delim('datasets/gsp1_emap_for_clustering.txt', delim = '\t', col_types = cols()) %>% 
#     gather(key = 'gene', value = score, -mutant) %>%
#     filter(gene %in% genes) %>%
#     spread(key = gene, value = score) %>% 
#     drop_na() # need to get rid of rows with all NA
  
#   mut.mat <- select(dset, -mutant) %>% as.matrix() %>% t()
#   colnames(mut.mat) <- dset %>% pull(mutant)
#   mut.hc <- clustfn(mut.mat)
#   arr.mat <- select(dset, -mutant) %>% as.matrix()
#   rownames(arr.mat) <- dset %>% pull(mutant)
#   arr.hc <- clustfn(arr.mat)
  
#   # GAP GEF dataset
#   dset_GG <- dset %>% filter(mutant %in% mutants_with_GAP_GEF_data)
#   mut_GG.mat <- select(dset_GG, -mutant) %>% as.matrix() %>% t()
#   colnames(mut_GG.mat) <- dset_GG %>% pull(mutant)
#   mut_GG.hc <- clustfn(mut_GG.mat)
#   arr_GG.mat <- select(dset_GG, -mutant) %>% as.matrix()
#   rownames(arr_GG.mat) <- dset_GG %>% pull(mutant)
#   arr_GG.hc <- clustfn(arr_GG.mat)

#   gl <- list()
  
#   gl[[1]] <- make_heatmap(mat = mut.mat.withNA, main = 'Full emap',
#                           Rowv = as.dendrogram(full_emap_rowtree),
#                           Colv = as.dendrogram(arr.hc))
#   gl[[2]] <- make_heatmap(mat = mut.mat.withNA, main = 'Seq order',
#                           Rowv = FALSE,
#                           Colv = as.dendrogram(arr.hc))
#   gl[[3]] <- make_heatmap(mat = arr.mat, main = 'Gene set recluster',
#                           Rowv = as.dendrogram(mut.hc),
#                           Colv = as.dendrogram(arr.hc))
#   gl[[4]] <- make_heatmap(mat = mut.mat.withNA, main = 'Interface',
#                           Rowv = as.dendrogram(interface_rowtree),
#                           Colv = as.dendrogram(arr.hc))
#   gl[[5]] <- grid.text(gene_set_list[i], gp=gpar(fontsize=26))
#   gl[[6]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Full emap',
#                           Rowv = as.dendrogram(full_emap_GAP_GEF_rowtree),
#                           Colv = as.dendrogram(arr_GG.hc))
#   gl[[7]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Seq order',
#                           Rowv = FALSE,
#                           Colv = as.dendrogram(arr_GG.hc))
#   gl[[8]] <- make_heatmap(mat = arr_GG.mat, main = 'Gene set recluster',
#                           Rowv = as.dendrogram(mut_GG.hc),
#                           Colv = as.dendrogram(arr_GG.hc))
#   gl[[9]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'Interface',
#                           Rowv = as.dendrogram(interface_GAP_GEF_rowtree),
#                           Colv = as.dendrogram(arr_GG.hc))
#   gl[[10]] <- make_heatmap(mat = mut_GG.mat.withNA, main = 'GAP GEF',
#                           Rowv = as.dendrogram(GAP_GEF_rowtree),
#                           Colv = as.dendrogram(arr_GG.hc))

#   grid.arrange(grobs=gl, ncol=5, clip=TRUE, newpage=TRUE)
#   plots_list[[i]] <- recordPlot()
# }
# pdf('images/heatmaps_corr.pdf', height = 20, width = 32, onefile = TRUE)
# for (my.plot in plots_list) {
#   replayPlot(my.plot)
# }
# graphics.off()