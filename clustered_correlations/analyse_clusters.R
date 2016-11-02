### This script can analyse the GO slim terms of the clusters made by hclust_EMAP_library_genes.R
## and also make it's own library clusters based on GO slim terms

library(tools)
library(magrittr)
library(tidyr)
options( stringsAsFactors = F)

setwd( "~/Documents/GSP1/E_MAP_data/emap_analysis/clustered_correlations/" )

clusters_prefix_to_analyse <- "2016-11-01"

GO_slims <- read.delim( "go_slim_mapping.tab.txt", head = F)
names(GO_slims) <- c("ORF", "Gene", "SGID", "GO_Aspect", "GO_Slim_term",
                     "GOID", "Feature_type")
GO_slims <- GO_slims[ GO_slims$GO_Aspect == "P",][,c(1,2,5)]

#### check for overly general slims
GO_slims_count <- with(GO_slims, aggregate(ORF, by = list(GO_Slim_term = GO_Slim_term), length ) )
names(GO_slims_count)[2] <- "GO_Slim_term_count"
GO_slims_count <- GO_slims_count[ order(GO_slims_count$GO_Slim_term_count, decreasing = T), ]
### remove some terms
GO_slim_terms_to_remove <- as.character(expression(
  signaling, biological_process, response, "response to chemical", "protein complex biogenesis", other
))
GO_slims <- GO_slims[ ! GO_slims$GO_Slim_term %in% GO_slim_terms_to_remove, ]

#########################
files_list <- list.files("clusters")
clusters_and_GO <- data.frame()
for ( f in seq_along(files_list) ) {
  file <- files_list[f]
  if ( grepl( clusters_prefix_to_analyse, file) & ! grepl(".pdf", file) & ! grepl("GO", file)) {
    
    cluster_method <- paste(unlist(strsplit(file_path_sans_ext(file), split = "_"))[2:4], collapse  = "_")
    
    cluster_df <- read.delim(file = file.path("clusters", file), head = T)
    if ( grepl("gene", file)) {
      cluster_df <- cluster_df[, c(1, 3, 4)]
    }
    GO_clusters_merged <- merge( cluster_df, GO_slims, by = "ORF")[, c(1, 2, 4, 5)]
    GO_clusters_merged <- GO_clusters_merged[order(GO_clusters_merged$cluster, GO_clusters_merged$Gene), ]
    by_cluster_GO_term_count <- with( GO_clusters_merged, aggregate( Gene, 
                by = list( cluster = cluster, GO_Slim_term = GO_Slim_term), length ) )
    names(by_cluster_GO_term_count)[3] <- "Slim_term_counts_per_cluster"
    cluster_sizes <- with( GO_clusters_merged, aggregate(ORF, by = list(cluster = cluster), length ) )
    names(cluster_sizes)[2] <- "cluster_size"
    by_cluster_GO_term_count <- merge(by_cluster_GO_term_count, cluster_sizes, by = "cluster")
    by_cluster_GO_term_count <- merge(by_cluster_GO_term_count, GO_slims_count, by = "GO_Slim_term")
    by_cluster_GO_term_count <- cbind( by_cluster_GO_term_count,
         data.frame( "percent_cluster" = by_cluster_GO_term_count$Slim_term_counts_per_cluster / by_cluster_GO_term_count$cluster_size, 
                     "percent_slim" = by_cluster_GO_term_count$Slim_term_counts_per_cluster / by_cluster_GO_term_count$GO_Slim_term_count,
                     "normalized" = by_cluster_GO_term_count$Slim_term_counts_per_cluster / (
                      by_cluster_GO_term_count$cluster_size * by_cluster_GO_term_count$GO_Slim_term_count
                      )
                    ))
    by_cluster_GO_term_count <- by_cluster_GO_term_count[ order( by_cluster_GO_term_count$cluster, 
                                      by_cluster_GO_term_count$normalized, decreasing = T ), ]

    filename = paste0(Sys.Date(), "_", file_path_sans_ext(file), "_", "_cluster_GO_terms.txt")
    write.table(by_cluster_GO_term_count, file = filename, quote = F, row.names = F, sep = "\t")
    
    by_cluster_GO_term_count <- cbind(by_cluster_GO_term_count, "clustering_method" = cluster_method)
    
    clusters_and_GO <- rbind(clusters_and_GO, by_cluster_GO_term_count)
  }
}

clusters_and_GO <- cbind( clusters_and_GO, "scaled_normalized" = scale (clusters_and_GO$normalized) )
clusters_and_GO <- clusters_and_GO[ order( clusters_and_GO$clustering_method, clusters_and_GO$cluster, decreasing = T), ]
plot(density(clusters_and_GO$scaled_normalized))
clusters_and_GO_filtered <- clusters_and_GO[ clusters_and_GO$Slim_term_counts_per_cluster > 2, ]
clusters_and_GO_filtered <- clusters_and_GO_filtered[ order( clusters_and_GO_filtered$clustering_method, clusters_and_GO_filtered$cluster), ]

# temp <- clusters_and_GO_filtered[clusters_and_GO_filtered$clustering_method == "library_pearson_complete",]
# write.table(temp, file = paste(Sys.Date(), "_gene_correlation_complete_cluster_GO_descriptions.txt"), quote = F, sep = "\t", row.names = F)



