library(gplots)

options(stringsAsFactors = F)

ms_partners_df <- read.delim("Gsp1_MS_hits.txt", head = F)
ms_partners <- ms_partners_df$V1
#gtpases_df <- read.delim("GTPase_list_uniq.txt", head = T)
#gtpases <- gtpases_df$V1
load("corr_of_corr.RData")
corr_of_corr_df <- result[["corr_of_corr"]]
clusters <- result[["clusters"]]
rm(result)
# corr_of_corr_df <- data.frame()
# clusters <- list()
# 
# dir_name <- "corr_of_corr"
# files_list <- list.files(dir_name)
# for ( i in seq_along(files_list) ) {
#   filename <- files_list[i]
#   if ( grepl("_corr_of_corr.RData", filename) ) {
#     file_path <- file.path(dir_name, filename)
#     load(file_path)
#     cluster <- final.list[["cluster"]]
#     clusters[[final.list[["cluster_number"]]]] <- cluster
#     corr_of_corr_df <- rbind(corr_of_corr_df, cbind(
#       final.list[["corr_of_corr_df"]], "cluster_number" = final.list[["cluster_number"]], "cluster" = cluster) )
#   }
# }
# 
# #check <- final.list[["corr_of_corr_df"]][order(final.list[["corr_of_corr_df"]]$mutant, 
#  #                           final.list[["corr_of_corr_df"]]$partner, final.list[["corr_of_corr_df"]]$partner),]
# #check <- head(check, n=1000)
# rm(final.list)
# 
# write.table(corr_of_corr_df, file = "GO_cluster_all_corr_of_corr.txt", quote = F, sep = "\t", row.names = F)

corr_of_corr_df <- data.frame(lapply(corr_of_corr_df, gsub, pattern = " - [A-z0-9]{7}", replacement = "", perl = T))
corr_of_corr_df <- corr_of_corr_df[ order(corr_of_corr_df$mutant, corr_of_corr_df$partner, corr_of_corr_df$method),]
partners.list<- list()
partners.list[["Partners_with_structures"]] <- c("RNA1", "SRM1", "YRB1", "PSE1", "MOG1")
partners.list[["MS_partners"]] <- ms_partners
partners.list[["All_Partners"]] <- as.character(unique(corr_of_corr_df[["partner"]]))
partners.list[["GTPases"]] <- gtpases
selected_measures <- c( "corr_of_corr_no_na_no_mut", "dot_product_no_na_no_mut")

temp <- corr_of_corr_df[ grepl(pattern = paste(unlist(partners.list[[1]]), collapse="|"), corr_of_corr_df[["partner"]]), ]
partners_to_color <- as.character(unique(temp[["partner"]]))
partners_to_color <- partners_to_color[! partners_to_color == "RNA15_TSQ652"]
coloring_list <- list()
colors <- rainbow(length(partners_to_color))
for ( i in seq_along(partners_to_color) ) {
  coloring_list[partners_to_color[i]] <- colors[i]
}



for ( m in seq_along(selected_measures) ) {
  measure = selected_measures[m]
  temp_measure_df <- corr_of_corr_df[
    corr_of_corr_df[["method"]] == measure, ]
  #merged_subclusters_corr_of_corr <- list()
  for ( partner in seq_along(ms_partners) ) {
    cluster <- clusters[[partner]]
    #print(cluster)
    #cluster_name_split <- unlist(strsplit(cluster, "_"))
    ##### group subclusters into a GO cluster, and group them all together
    temp.cluster.df <- temp_measure_df[ temp_measure_df[["cluster_number"]] == partner, ]
    temp.cluster.df <- temp.cluster.df[ order(temp.cluster.df$mutant, temp.cluster.df$partner), ]
    #merged_subclusters_corr_of_corr <- rbind(merged_subclusters_corr_of_corr,
    #       cbind(temp.cluster.df[,1:4], "GO_cluster" = cluster_name_split[1], "subcluster" = cluster_name_split[3]))
    
    for ( l in seq_along( partners.list ) ) {
      
      temp.df <- temp.cluster.df[ grepl(pattern = paste(unlist(partners.list[[l]]), collapse="|"), temp.cluster.df[["partner"]]), ]
      #temp.df <- temp.cluster.df[ temp.cluster.df[["partner"]] %in% unlist( partners.list[[l]] ), ]  
      
      if (l == 1) {
        temp.df <- temp.df[! temp.df$partner == "RNA15_TSQ652",]
      }
      partners.to.remove <- unique(temp.df$partner[temp.df$value == "NA"])
      if (length(partners.to.remove) > 0) {
        temp.df <- temp.df[! grepl(pattern = paste(partners.to.remove, collapse="|"), temp.df[["partner"]]), ]
      }
      temp.partners <- as.character( unique( temp.df$partner ))
      corr.mutants <- as.character(unique( temp.df[["mutant"]] ))
      temp.df <- temp.df[ order(temp.df$mutant, temp.df$partner), ]
      if ( nrow(temp.df) > 2) {
        correlation.matrix <- matrix( as.numeric(temp.df$value), byrow = T, nrow = length(corr.mutants), ncol = length(temp.partners) )
        correlation.matrix <- t(correlation.matrix)
        colnames(correlation.matrix) <- corr.mutants
        rownames(correlation.matrix) <- as.character( unique( temp.df$partner ) ) 
        filename <- paste0(Sys.Date(), "_", names(partners.list[l]), "_", measure, "_", cluster, ".pdf")
        print(measure)
        print(cluster)
        print(names(partners.list[l]))
        pdf(file = filename, width = 14, height = 14)
        
        cols <- rep('black', nrow(correlation.matrix))
        for (i in seq_along(names(coloring_list))) {
          cols[row.names(correlation.matrix) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
        }
        heatmap <- heatmap.2( correlation.matrix, scale = "none", dendrogram="both", 
                              trace="none", density.info="none", col = cm.colors, 
                              key.ylab="partners", key.xlab = m, 
                              main = names(partners.list[l]), cexCol=0.9, cexRow = 0.9, keysize=0.9, margin=c(10,10), colRow = cols)
        
        dev.off()
      }
    }
  }
}






# for ( m in seq_along(selected_measures) ) {
#   measure = selected_measures[m]
#   temp_measure_df <- corr_of_corr_df[
#     corr_of_corr_df[["method"]] == measure, ]
#   #merged_subclusters_corr_of_corr <- list()
#   for ( clust_n in seq_along(clusters) ) {
#     cluster <- clusters[[clust_n]]
#     #print(cluster)
#     #cluster_name_split <- unlist(strsplit(cluster, "_"))
#     ##### group subclusters into a GO cluster, and group them all together
#     temp.cluster.df <- temp_measure_df[ temp_measure_df[["cluster_number"]] == clust_n, ]
#     temp.cluster.df <- temp.cluster.df[ order(temp.cluster.df$mutant, temp.cluster.df$partner), ]
#     #merged_subclusters_corr_of_corr <- rbind(merged_subclusters_corr_of_corr,
#      #       cbind(temp.cluster.df[,1:4], "GO_cluster" = cluster_name_split[1], "subcluster" = cluster_name_split[3]))
#     
#     for ( l in seq_along( partners.list ) ) {
#       
#       temp.df <- temp.cluster.df[ grepl(pattern = paste(unlist(partners.list[[l]]), collapse="|"), temp.cluster.df[["partner"]]), ]
#       #temp.df <- temp.cluster.df[ temp.cluster.df[["partner"]] %in% unlist( partners.list[[l]] ), ]  
#       
#       if (l == 1) {
#         temp.df <- temp.df[! temp.df$partner == "RNA15_TSQ652",]
#       }
#       partners.to.remove <- unique(temp.df$partner[temp.df$value == "NA"])
#       if (length(partners.to.remove) > 0) {
#           temp.df <- temp.df[! grepl(pattern = paste(partners.to.remove, collapse="|"), temp.df[["partner"]]), ]
#       }
#       temp.partners <- as.character( unique( temp.df$partner ))
#       corr.mutants <- as.character(unique( temp.df[["mutant"]] ))
#       temp.df <- temp.df[ order(temp.df$mutant, temp.df$partner), ]
#       if ( nrow(temp.df) > 2) {
#         correlation.matrix <- matrix( as.numeric(temp.df$value), byrow = T, nrow = length(corr.mutants), ncol = length(temp.partners) )
#         correlation.matrix <- t(correlation.matrix)
#         colnames(correlation.matrix) <- corr.mutants
#         rownames(correlation.matrix) <- as.character( unique( temp.df$partner ) ) 
#         filename <- paste0(Sys.Date(), "_", names(partners.list[l]), "_", measure, "_", cluster, ".pdf")
#         print(measure)
#         print(cluster)
#         print(names(partners.list[l]))
#         pdf(file = filename, width = 14, height = 14)
#         
#         cols <- rep('black', nrow(correlation.matrix))
#         for (i in seq_along(names(coloring_list))) {
#           cols[row.names(correlation.matrix) %in% as.vector(names(coloring_list)[i]) ] <- unlist(coloring_list[[i]])
#         }
#         heatmap <- heatmap.2( correlation.matrix, scale = "none", dendrogram="both", 
#                               trace="none", density.info="none", col = cm.colors, 
#                               key.ylab="partners", key.xlab = m, 
#                               main = names(partners.list[l]), cexCol=0.9, cexRow = 0.9, keysize=0.9, margin=c(10,10), colRow = cols)
#         
#         dev.off()
#       }
#     }
#   }
# }
# 
