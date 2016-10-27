library(fastcluster)  # overwrites the hclust function and makes it faster
library(dynamicTreeCut)
library(amap)

setwd("~/Documents/GSP1/E_MAP_data/June2016_analysis/clustered_correlations")

options(stringsAsFactors = FALSE)

#### load an orf to gene name index
orf_gene_name_index <- read.delim("orf_gene_GO_sgd_annotation.txt", head = F)
orf_index <- unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)

###### inputting data
ubermap <- read.delim("preprocessed_ubermap_ubergenes_only.txt", head = T)  # preprocessed file made by preprocess_ubermap_merged_data.R
#ubermap <- read.delim(opt$ubermap, head = T)
#head(ubermap)
#  library.ORF Gene.ORF     score Gene_uniq Gene.gene_name library.gene_name
### this is the test code to check that everything is nr
# count.control.Gene <- with(ubermap, aggregate(Gene_uniq, by = list(library.ORF = library.ORF), length))
# range(count.control.Gene$x)
# count.control.library <- with(ubermap, aggregate(library.ORF, by = list(Gene_uniq = Gene_uniq), length))
# count.control.library <- count.control.library[order(count.control.library$x, decreasing = T),]
# unique(count.control.library$Gene_uniq[count.control.library$x > 1354])
# temp.ubermap <- subset(ubermap, Gene_uniq == "YDR113C")
# range(count.control.library$x)         
mut.emap <- read.delim("preproceessed_ubermap_mut_only.txt", head = T)
#mut.emap <- read.delim(opt$emap, head = T)
### the ubermap is ordered first by Gene_uniq and then by library.ORF


#### cluster the library genes
### by pearson correlation, and by e-map score
### make a matrix
lib.ubermap.matrix <- list()  # all the matrices are stored in a single list

### library genes, e-map score only matrix
Gene_uniq <- as.character(unique(ubermap$Gene_uniq))
library.ORF <- as.character(unique(ubermap$library.ORF))
lib.ubermap.matrix[["library"]] <- matrix(ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(Gene_uniq))
rownames(lib.ubermap.matrix[["library"]]) <- library.ORF
colnames(lib.ubermap.matrix[["library"]]) <- Gene_uniq
#head(lib.ubermap.matrix[["library"]])[,1:20]  ### rows are library gene ORFs and columns are unique Gene identifiers (check preprocess_ubermap_merged_data.R for why uniq)

# library genes, pearson correlations matrix
# lib.ubermap.matrix[["library_pearson"]] <- cor(t(lib.ubermap.matrix[["library_e-map"]]), method = "pearson", use="pairwise.complete.obs")
# rownames(lib.ubermap.matrix[["library_pearson"]]) <- library.ORF
# colnames(lib.ubermap.matrix[["library_pearson"]]) <- library.ORF

#### Also do the clustering by Genes (this could be used to filter the library clustering further,
### so only consider clustering those query genes that also exist as library genes)
####   this ubermap subset only keeps the rows where the Gene.ORF also exists as a library ORF 
GeneORF_subset_ubermap <- ubermap[ubermap$Gene.ORF %in% library.ORF,]
#GeneORF_ubermap_library_ORFs <- as.character(unique(GeneORF_subset_ubermap$library.ORF))  ### this is the same as library.ORF
GeneORF_Gene_uniq <- as.character(unique(GeneORF_subset_ubermap$Gene_uniq))
Gene_uniq_Gene_ORF_pairs <- data.frame("Gene.ORF" = GeneORF_subset_ubermap$Gene.ORF, "Gene_uniq" = GeneORF_subset_ubermap$Gene_uniq)
Gene_uniq_Gene_ORF_pairs <- Gene_uniq_Gene_ORF_pairs[! duplicated(Gene_uniq_Gene_ORF_pairs),]
### make the matrix as above and then trasnpose it - it's faster than reordering the data.frame
#### query genes, e-map scores
lib.ubermap.matrix[["gene"]] <- matrix(GeneORF_subset_ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(GeneORF_Gene_uniq))
#head(lib.ubermap.matrix[["gene"]])[,1:20]
rownames(lib.ubermap.matrix[["gene"]]) <- library.ORF
colnames(lib.ubermap.matrix[["gene"]]) <- GeneORF_Gene_uniq
lib.ubermap.matrix[["gene"]] <- t(lib.ubermap.matrix[["gene"]])

# query genes, pearson correlations
# lib.ubermap.matrix[["gene_pearson"]] <- cor(t(lib.ubermap.matrix[["gene"]]), method = "pearson", use = "pairwise.complete.obs")
# rownames(lib.ubermap.matrix[["gene_pearson"]]) <- GeneORF_Gene_uniq
# colnames(lib.ubermap.matrix[["gene_pearson"]]) <- GeneORF_Gene_uniq


distance_methods <- c("pearson", "correlation", "euclidean")
hclust_methods <- c("complete", "average")
### do the clustering for each of the matrices
clustering_summary_df <- data.frame()
pdf(paste0(Sys.Date(), "_clusters.pdf"), width = 10)
for (mat in names(lib.ubermap.matrix)) {
  temp.mat <- lib.ubermap.matrix[[mat]]
  for (i in seq_along(distance_methods)) {
    for (j in seq_along(hclust_methods)) {
      dist_method <- distance_methods[i]
      hclust_method <- hclust_methods[j]
      dissim <- Dist(temp.mat, method = dist_method)
      dendro <- hclust(dissim, method = hclust_method)
      plot(dendro, cex = 0.2, main = paste(mat, dist_method, hclust_method, sep = "_"))
      cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, minClusterSize = 30, method = "hybrid", deepSplit = 4, pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
      clusters <- data.frame(cut_hybrid, rownames(temp.mat))
      if (grepl("gene", mat)) {
        names(clusters) <- c('cluster', 'Gene_uniq')
        clusters <- merge(clusters, Gene_uniq_Gene_ORF_pairs, by = "Gene_uniq")
        clusters <- merge(clusters, orf_index, by.x = "Gene.ORF", by.y = "orf")
        names(clusters)[1] <- "ORF"
        clusters <- clusters[order(clusters$cluster, decreasing = T),]
      } else {
        names(clusters) <- c('cluster', 'ORF')
        clusters <- merge(clusters, orf_index, by.x = "ORF", by.y = "orf")
        clusters <- clusters[order(clusters$cluster, decreasing = T),]
      }
      clustering_summary_df <- rbind(clustering_summary_df,
                                     data.frame("data" = mat, "dist_method" = dist_method, "hclust_method" = hclust_method,
                                     "n_clusters" = max(clusters$cluster), 
                                     "number_of_unclustered_genes" = length(clusters$ORF[clusters$cluster == 0])))
      outfilename <- paste(Sys.Date(), mat, dist_method, hclust_method, "clusters.txt", sep = "_")
      write.table(clusters, file = outfilename, quote = F, row.names = F, sep = "\t")
    }
  }
}
dev.off()


# ##### cluster the Gene_uniq (and process it afterwards to get ORFs for merging with library.clusters)
# Gene_dissim <- dist(GeneORF.ubermap.matrix, method = "euclidean")
# Gene_dendro <- hclust(Gene_dissim, method = "complete")
# plot(Gene_dendro, cex = 0.2)
# Gene_cut_hybrid <- cutreeDynamic(dendro = Gene_dendro, cutHeight = NULL, minClusterSize = 30, method = "hybrid", deepSplit = 4, pamStage = T, distM = as.matrix(Gene_dissim), maxPamDist = 0, verbose = 0)
# Gene.clusters <- data.frame(Gene_cut_hybrid, rownames(GeneORF.ubermap.matrix))
# Gene.clusters <- Gene.clusters[order(Gene.clusters$Gene_cut_hybrid, decreasing = T),]
# head(Gene.clusters)  ## 19 clusters
# length(Gene.clusters$Gene_cut_hybrid[Gene.clusters$Gene_cut_hybrid == 0]) ### 222 genes not in any cluster
# names(Gene.clusters) <- c('cluster', 'Gene_uniq')
# Gene_uniq_Gene_ORF_pairs <- data.frame("Gene.ORF" = GeneORF_subset_ubermap$Gene.ORF, "Gene_uniq" = GeneORF_subset_ubermap$Gene_uniq)
# Gene_uniq_Gene_ORF_pairs <- Gene_uniq_Gene_ORF_pairs[! duplicated(Gene_uniq_Gene_ORF_pairs),]
# Gene.clusters <- merge(Gene.clusters, Gene_uniq_Gene_ORF_pairs, by = "Gene_uniq")
# Gene.clusters <- merge(Gene.clusters, orf_index, by.x = "Gene.ORF", by.y = "orf")
# Gene.clusters <- Gene.clusters[order(Gene.clusters$cluster, decreasing = T),]
# 
# ### clusters are named as numbers - get overlaps between library and Gene clusters to get matching cluster names (numbers)
# cluster.overlaps <- data.frame()
# for (gc in 1:max(Gene.clusters$cluster)) {
#   for (lc in 1:max(library.clusters$cluster)) {
#     temp.gene.cluster <- Gene.clusters$Gene.ORF[Gene.clusters$cluster == gc]
#     temp.library.cluster <- library.clusters$library_gene_ORF[library.clusters$cluster == lc]
#     overlap <- get_overlap(temp.gene.cluster, temp.library.cluster)
#     cluster.overlaps <- rbind(cluster.overlaps, 
#         data.frame("gene_cluster" = gc, "gene_cluster_size" = length(temp.gene.cluster),
#                    "library_cluster" = lc, "library_cluster_size" = length(temp.library.cluster),
#                    "overlap" = overlap
#                   ))
#   }
# }
# cluster.overlaps <- cluster.overlaps[order(cluster.overlaps$gene_cluster, cluster.overlaps$overlap, decreasing = T),]
# 
# 
# ### for each pair of genes in the same cluster in Gene.clusters check if they are also in the same cluster in library.clusters
# ### if yes, keep that pair in that cluster
# 
# 
# write.table(library.clusters, file = "library_genes_hclusted.txt", quote = F, sep = "\t", row.names = F)
# 
# 
