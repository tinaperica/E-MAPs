options(stringsAsFactors = F)
# run a task per cluster - there are 35 clusters for the GO_clustering approach (and also corresponding 35 random clusters)
cluster_number <- as.numeric( Sys.getenv( "SGE_TASK_ID" ))
#cluster_number <- 1
#### preprocessed correlations
load( paste0( "preprocessed_correlations_random/", "random_clusters_corr_", cluster_number, ".RData" ) )
outputfilename <- paste0("Output/random_clusters/random_cluster_", cluster_number, "_corr_of_corr.RData")

cluster <- corr_lst[["cluster"]]
mutants <- as.character(unique(corr_lst[["table_of_pairs"]][["mutant"]]))
###
mutant_partner_corr_of_corr <- data.frame()
for (p in seq_len( nrow(corr_lst[["table_of_pairs"]]) ) ) {
  mutant <- as.character( corr_lst[["table_of_pairs"]][["mutant"]][p] )
  partner <- as.character( corr_lst[["table_of_pairs"]][["partner"]][p] ) 
  mut_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene1"]] == mutant, ]
  temp_mut_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene2"]] == mutant, ]
  names(temp_mut_correlations_df) <- c("gene2", "gene1","cluster_method", "cluster", "corr")
  temp_mut_correlations_df <- temp_mut_correlations_df[,c("gene1", "gene2", "cluster_method", "cluster", "corr")]
  temp_mut_correlations_df <- temp_mut_correlations_df[, c(2,1,3:5)]
  mut_correlations_df <- rbind( mut_correlations_df, temp_mut_correlations_df)
  partner_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene1"]] == partner, ]
  temp_partner_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene2"]] == partner, ]
  names(temp_partner_correlations_df) <- c("gene2", "gene1", "cluster_method", "cluster", "corr")
  temp_partner_correlations_df <- temp_partner_correlations_df[, c(2, 1, 3:5)]
  partner_correlations_df <- rbind( partner_correlations_df, temp_partner_correlations_df)
  if ( (nrow(mut_correlations_df) > 3 ) & ( nrow(partner_correlations_df) > 3) ) { 
    merge_for_corr <- merge(mut_correlations_df, partner_correlations_df, by = c("gene2", "cluster", "cluster_method"))
    merge_for_corr <- merge_for_corr[! duplicated(merge_for_corr), ]
    merge_for_corr <- merge_for_corr[complete.cases(merge_for_corr),]
    temp.corr_of_corr <- data.frame()
    if ( nrow(merge_for_corr) > 3 ) {
      names(merge_for_corr) <- expression(uber_gene, cluster, cluster_method, mutant, mutant.corr, partner, partner.corr)
      #### calculate dot product and corr_of_corr with a series of filtering variations
      ### remove correlations with the mutants (remove if other_gene is a mutant) and repeat everything
      merge_for_corr_no_mut <- merge_for_corr[ ! merge_for_corr[["uber_gene"]] %in% mutants, ]
      if ( nrow(merge_for_corr_no_mut) > 3 ) {
        
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = 
                                                                   cor(merge_for_corr_no_mut$mutant.corr, merge_for_corr_no_mut$partner.corr)))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = 
                                                                   (merge_for_corr_no_mut$mutant.corr %*% merge_for_corr_no_mut$partner.corr) / length(merge_for_corr_no_mut$mutant.corr)))
        
      } else {
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))

      }
    } else {
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))
    }
  } else {
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))
  }
  temp.corr_of_corr <- cbind(temp.corr_of_corr, data.frame("mutant" = mutant, "partner" = partner))
  mutant_partner_corr_of_corr <- rbind(mutant_partner_corr_of_corr, temp.corr_of_corr)
}


final.list<-list()
final.list[["cluster"]] <- cluster
final.list[["cluster_number"]] <- cluster_number
final.list[["corr_of_corr_df"]] <- mutant_partner_corr_of_corr

save(final.list, file = outputfilename)
