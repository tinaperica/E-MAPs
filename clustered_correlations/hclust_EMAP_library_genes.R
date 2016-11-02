### this script makes hierarchical clusters of genes based on ubermap e-map scores
### as well as based on GO slims terms

library(fastcluster)  # overwrites the hclust function and makes it faster
library(dynamicTreeCut)
library(amap)
library(magrittr)

setwd("~/Documents/GSP1/E_MAP_data/emap_analysis/clustered_correlations")

options(stringsAsFactors = FALSE)

#### load an orf to gene name index
orf_gene_name_index <- read.delim("../orf_gene_GO_sgd_annotation.txt", head = F)
orf_index <- unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)

###### inputting ubermap data
#ubermap <- read.delim("preprocessed_ubermap_ubergenes_only.txt", head = T)  # preprocessed file made by preprocess_ubermap_merged_data.R
ubermap <- read.delim("../preprocessed_ubermap_ubergenes_only_significant.txt", head = T)
mut.emap <- read.delim("../preprocessed_ubermap_mut_only_significant.txt", head = T)

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
### the ubermap is ordered first by Gene_uniq and then by library.ORF


### also make clusters based on GO slims
GO_slims <- read.delim( "go_slim_mapping.tab.txt", head = F)
names(GO_slims) <- c("ORF", "Gene", "SGID", "GO_Aspect", "GO_Slim_term",
                     "GOID", "Feature_type")
#GO_slims <- GO_slims[ GO_slims$GO_Aspect == "P",][,c(1,2,5)]
GO_slims <- GO_slims[,c(1,2,5)]
#### check for overly general slims
GO_slims_count <- with(GO_slims, aggregate(ORF, by = list(GO_Slim_term = GO_Slim_term), length ) )
names(GO_slims_count)[2] <- "GO_Slim_term_count"
GO_slims_count <- GO_slims_count[ order(GO_slims_count$GO_Slim_term_count, decreasing = T), ]
### remove some terms
GO_slim_terms_to_remove <- as.character(expression(
  signaling, biological_process, response, "response to chemical", "protein complex biogenesis", other, "not_yet_annotated",
  "ion binding", "cellular_component", "molecular_function", "cytoplasm", "structural molecule activity", "ATPase activity"
))
GO_slims <- GO_slims[ ! GO_slims$GO_Slim_term %in% GO_slim_terms_to_remove, ]
##############################################
### manually select some terms based on which I will make GO based clusters - > this approach gives a small number of relativelly general clusters
GO_slims_for_clusters <- list()
GO_slims_for_clusters[["ribosome"]] <- unique ( GO_slims$GO_Slim_term[ 
  grep( paste( c("rRNA", "ribosom", "translation", "tRNA", "snoRNA"), collapse = "|"), GO_slims$GO_Slim_term)
  ])
GO_slims_for_clusters[["transcription and mRNA processing"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mRNA", "splicing", "transcription", "RNA modification"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["Golgi and ER"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("protein lipidation", "protein maturation", "endocytosis", "regulation of transport", "glycosylation",
                 "vesicle organization","endosom", "Golgi"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["peroxisome"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("peroxisom"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["vacuole"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("vacuol"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["mitochondrion"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mitochond", "respirat"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["chromatin"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("chromatin","histone", "telomere", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["cytoskeleton"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cytoskelet", "cytokinesis", "chromosome segregation", "microtubule"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["cell cycle"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cell cycle", "recombination", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["budding"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("bud", "fission", "sporulation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["lipids"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("lipid"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["nuclear transport"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("nucleus", "nuclear"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["metabolic"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("metabolic", "metabolite"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
#####################################################################
#### cluster the library genes
### by pearson correlation, and by e-map score
### also by GO slims cat
### make a matrix
lib.ubermap.matrix <- list()  # all the matrices are stored in a single list


#### cluster output files have format ORF cluster gene_name
broad_GO_cat_clusters <- data.frame()
for (cat in names(GO_slims_for_clusters)) {
  GO_slims_cat <- GO_slims_for_clusters[[cat]]
  temp_GO_slims <- GO_slims[ GO_slims$GO_Slim_term %in% GO_slims_cat, ]
  temp_GO_slim_ORFS <- as.character(unique(temp_GO_slims$ORF))
  temp_ubermap <- ubermap [ubermap$library.ORF %in% temp_GO_slim_ORFS,]
  clusters <- temp_ubermap[,c(1,6)]
  clusters <- clusters[ ! duplicated(clusters), ]
  clusters <- cbind(clusters, "cluster" = cat)
  names(clusters) <- c("ORF", "gene_name", "cluster")
  clusters <- clusters[, c(1,3,2)]
  broad_GO_cat_clusters <- rbind(broad_GO_cat_clusters, clusters)
  ### add each of these clusters onto the lib.ubermap.matrix so I can subscluster them with hclust later
  Gene_uniq <- as.character(unique(temp_ubermap$Gene_uniq))
  library.ORF <- as.character(unique(temp_ubermap$library.ORF))
  lib.ubermap.matrix[[paste0(cat, "_GO")]] <- matrix( temp_ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(Gene_uniq) )
  rownames( lib.ubermap.matrix[[ paste0( cat, "_GO" ) ]] ) <- library.ORF
  colnames( lib.ubermap.matrix[[ paste0( cat, "_GO" ) ]] ) <- Gene_uniq
  
}
outfilename <- paste(Sys.Date(), "broad_GO_cat_clusters.txt", sep = "_")
write.table(broad_GO_cat_clusters, file = outfilename, quote = F, row.names = F, sep = "\t")
head(lib.ubermap.matrix[["ribosome_GO"]])[,1:20]

# Also make clusters based on individual slim GO categories (that will be 100 clusters, but they will be highly overlapping)
## this is either a stupid idea or it will work best




### library genes, e-map score only matrix
Gene_uniq <- as.character(unique(ubermap$Gene_uniq))
library.ORF <- as.character(unique(ubermap$library.ORF))
lib.ubermap.matrix[["library"]] <- matrix(ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(Gene_uniq))
rownames(lib.ubermap.matrix[["library"]]) <- library.ORF
colnames(lib.ubermap.matrix[["library"]]) <- Gene_uniq
#head(lib.ubermap.matrix[["library"]])[,1:20]  ### rows are library gene ORFs and columns are unique Gene identifiers (check preprocess_ubermap_merged_data.R for why uniq)


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


distance_methods <- c("pearson", "correlation", "euclidean")
hclust_methods <- c("complete", "average")
### do the clustering for each of the matrices
clustering_summary_df <- data.frame()
broad_GO_cat_subclusters <- data.frame()
pdf(paste0(Sys.Date(), "_clusters.pdf"), width = 10)
for (i in seq_along(distance_methods)) {
  for (j in seq_along(hclust_methods)) {
    GO_slims_subcluster_unclustered_genes <- vector()
    GO_slims_subcluster_unique_clustered_genes <- vector()
    for (mat in names(lib.ubermap.matrix)) {
      temp.mat <- lib.ubermap.matrix[[mat]]
      dist_method <- distance_methods[i]
      hclust_method <- hclust_methods[j]
      dissim <- Dist(temp.mat, method = dist_method)
      dendro <- hclust(dissim, method = hclust_method)
      if (grepl("GO", mat)) {
        plot(dendro, cex = 0.2, main = paste(mat, dist_method, hclust_method, sep = "_"))
        cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, minClusterSize = 30, method = "hybrid", deepSplit = 4, pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
        clusters <- data.frame(cut_hybrid, rownames(temp.mat))
        names(clusters) <- c('subcluster', 'ORF')
        clusters <- cbind( clusters, "GO" = mat)
        if ( length( clusters[ clusters$subcluster != 0, ][,1] ) > 0 ) {  ## if there are not enough genes to make multiple clusters of min size 15, they will all be in one unclustered cluster 0
          ### in that case just treat that single 0 cluster as a cluster
          ### but if there are normal clusters, remove the 0 cluster
          GO_slims_subcluster_unclustered_genes <- append(GO_slims_subcluster_unclustered_genes, clusters$ORF[clusters$subcluster == 0])
          clusters <- clusters[ ! clusters$subcluster == 0, ]
        }
        GO_slims_subcluster_unique_clustered_genes <- append(GO_slims_subcluster_unique_clustered_genes, unique(clusters$ORF))
        clusters <- cbind( clusters, "cluster" = paste(clusters$GO, clusters$subcluster, sep = "_"))
        clusters <- clusters[,c(4,2)]
        clusters <- merge(clusters, orf_index, by.x = "ORF", by.y = "orf")
        clusters <- clusters[order(clusters$cluster, decreasing = T),]
        broad_GO_cat_subclusters <- rbind(broad_GO_cat_subclusters, clusters)
      } else {
        plot(dendro, cex = 0.2, main = paste(mat, dist_method, hclust_method, sep = "_"))
        cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, minClusterSize = 30, method = "hybrid", deepSplit = 4, pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
        clusters <- data.frame(cut_hybrid, rownames(temp.mat))
        if (grepl("gene", mat)) {
          names(clusters) <- c('cluster', 'Gene_uniq')
          clusters <- merge(clusters, Gene_uniq_Gene_ORF_pairs, by = "Gene_uniq")
          clusters <- merge(clusters, orf_index, by.x = "Gene.ORF", by.y = "orf")
          names(clusters)[1] <- "ORF"
          clusters <- clusters[order(clusters$cluster, decreasing = T),]
        } else if (grepl("library", mat)) {
          names(clusters) <- c('cluster', 'ORF')
          clusters <- merge(clusters, orf_index, by.x = "ORF", by.y = "orf")
          clusters <- clusters[order(clusters$cluster, decreasing = T),]
        } 
        clustering_summary_df <- rbind(clustering_summary_df,
                                       data.frame("data" = mat, "dist_method" = dist_method, "hclust_method" = hclust_method,
                                                  "n_clusters" = max(clusters$cluster), 
                                                  "number_of_unclustered_genes" = length(clusters$ORF[clusters$cluster == 0]),
                                                  "number_of_unique_clustered_genes" = length(clusters$ORF[clusters$cluster != 0])
                                                  ))
        outfilename <- paste(Sys.Date(), mat, dist_method, hclust_method, "clusters.txt", sep = "_")
        write.table(clusters, file = outfilename, quote = F, row.names = F, sep = "\t")
      }
    }
    clustering_summary_df <- rbind(clustering_summary_df,
                                   data.frame("data" = "GO_slims", "dist_method" = dist_method, "hclust_method" = hclust_method,
                                              "n_clusters" = length(unique(broad_GO_cat_subclusters$cluster)), 
                                              "number_of_unclustered_genes" = length( unique(GO_slims_subcluster_unclustered_genes)),
                                              "number_of_unique_clustered_genes" = length( unique (GO_slims_subcluster_unique_clustered_genes))
                                              ))
    
    outfilename <- paste(Sys.Date(), "GO_slims", dist_method, hclust_method, "clusters.txt", sep = "_")
    write.table(broad_GO_cat_subclusters, file = outfilename, quote = F, row.names = F, sep = "\t")
  }
}
dev.off()




