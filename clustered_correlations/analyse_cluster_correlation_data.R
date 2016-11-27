options(stringsAsFactors = F)
# run a task per cluster - there are 44 clusters for the GO_clustering approach
cluster_number <- as.numeric( Sys.getenv( "SGE_TASK_ID" ))
#cluster_number <- 1
#### preprocessed correlations
load( paste0( "preprocessed_correlations/", "2016-11-26_", cluster_number, "_corr_lst.RData" ) )
cluster <- corr_lst[["cluster"]]
mutants <- as.character(unique(corr_lst[["table_of_pairs"]][["mutant"]]))
### get only part of the 
mutant_partner_corr_of_corr <- data.frame()
for (p in seq_len( nrow(corr_lst[["table_of_pairs"]]) ) ) {
  mutant <- as.character( corr_lst[["table_of_pairs"]][["mutant"]][p] )
  partner <- as.character( corr_lst[["table_of_pairs"]][["partner"]][p] ) 
  mut_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene1"]] == mutant, ]
  temp_mut_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene2"]] == mutant, ]
  names(temp_mut_correlations_df) <- c("gene2", "gene1", "cluster", "corr", "filter")
  temp_mut_correlations_df <- temp_mut_correlations_df[,c("gene1", "gene2", "cluster", "corr", "filter")]
  mut_correlations_df <- rbind( mut_correlations_df, temp_mut_correlations_df)
  partner_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene1"]] == partner, ]
  temp_partner_correlations_df <- corr_lst[["correlations"]][
    corr_lst[["correlations"]][["gene2"]] == partner, ]
  names(temp_partner_correlations_df) <- c("gene2", "gene1", "cluster", "corr", "filter")
  partner_correlations_df <- rbind( partner_correlations_df, temp_partner_correlations_df)
  if ( (nrow(mut_correlations_df) != 0 ) & ( nrow(partner_correlations_df) != 0) ) { 
    merge_for_corr <- merge(mut_correlations_df, partner_correlations_df, by = c("gene2", "cluster"))
    merge_for_corr <- merge_for_corr[! duplicated(merge_for_corr), ]
    merge_for_corr <- merge_for_corr[complete.cases(merge_for_corr),]
    temp.corr_of_corr <- data.frame()
    if ( nrow(merge_for_corr) > 3 ) {
      names(merge_for_corr) <- expression(other_gene, cluster, mutant, mutant.corr, filter.mut, partner, partner.corr, filter.partner)
      #### calculate dot product and corr_of_corr with a series of filtering variations
      #### treating all correlations as valid
      
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_all", "value" = 
        cor(merge_for_corr$mutant.corr, merge_for_corr$partner.corr)))
      
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_all", "value" = 
      (merge_for_corr$mutant.corr %*% merge_for_corr$partner.corr) / length(merge_for_corr$mutant.corr)))
      
      #### remove NA filter
      no_na_for_corr <- merge_for_corr[! (merge_for_corr$filter.mut == "NA" | merge_for_corr$filter.partner == "NA"), ]
      if (nrow(no_na_for_corr) > 3) {
        
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na", "value" =
          cor(no_na_for_corr$mutant.corr, no_na_for_corr$partner.corr)))
        
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na", "value" =
          (no_na_for_corr$mutant.corr %*% no_na_for_corr$partner.corr) / length(no_na_for_corr$mutant.corr)))
        
        ### remove all unfiltered (as well as no NA)
        filtered_for_corr <- no_na_for_corr[no_na_for_corr$filter.mut == "filtered" & no_na_for_corr$filter.partner == "filtered",]
        if ( nrow(filtered_for_corr) > 3 ) {
          
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered", "value" =
                                                                     cor(filtered_for_corr$mutant.corr, filtered_for_corr$partner.corr)))
          
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered", "value" =
                                                                     (filtered_for_corr$mutant.corr %*% filtered_for_corr$partner.corr) / length(filtered_for_corr$mutant.corr)))
          
        } else {
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered", "value" = "NA"))
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered", "value" = "NA"))
        }
      
      } else {
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered", "value" = "NA"))
      }
      
      ### remove correlations with the mutants (remove if other_gene is a mutant) and repeat everything
      merge_for_corr_no_mut <- merge_for_corr[ ! merge_for_corr[["other_gene"]] %in% mutants, ]
      if ( nrow(merge_for_corr_no_mut) > 5 ) {
        
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = 
              cor(merge_for_corr_no_mut$mutant.corr, merge_for_corr_no_mut$partner.corr)))
        
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = 
              (merge_for_corr_no_mut$mutant.corr %*% merge_for_corr_no_mut$partner.corr) / length(merge_for_corr_no_mut$mutant.corr)))
        
        no_na_for_corr_no_mut <- merge_for_corr_no_mut[! (merge_for_corr_no_mut$filter.mut == "NA" | merge_for_corr_no_mut$filter.partner == "NA"), ]
        if ( nrow(no_na_for_corr_no_mut) > 3 ) {
          
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na_no_mut", "value" =
                              cor(no_na_for_corr_no_mut$mutant.corr, no_na_for_corr_no_mut$partner.corr)))
          
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na_no_mut", "value" =
                        (no_na_for_corr_no_mut$mutant.corr %*% no_na_for_corr_no_mut$partner.corr) / length(no_na_for_corr_no_mut$mutant.corr)))
          
          ### remove all unfiltered (as well as no NA)
          filtered_for_corr_no_mut <- no_na_for_corr_no_mut[no_na_for_corr_no_mut$filter.mut == "filtered" & no_na_for_corr_no_mut$filter.partner == "filtered",]
          if ( nrow(filtered_for_corr_no_mut) > 3 ) {
            
            temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" =
                                cor(filtered_for_corr_no_mut$mutant.corr, filtered_for_corr_no_mut$partner.corr)))
            
            temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" =
                                (filtered_for_corr_no_mut$mutant.corr %*% filtered_for_corr_no_mut$partner.corr) / length(filtered_for_corr_no_mut$mutant.corr)))
            
          } else {
            temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" = "NA"))
            temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" = "NA"))
          }
          
        } else {
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na_no_mut", "value" = "NA"))
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na_no_mut", "value" = "NA"))
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" = "NA"))
          temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" = "NA"))
        }

      } else {
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" = "NA"))
        temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" = "NA"))
        
      }
      
    } else {
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_all", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_all", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered", "value" = "NA"))
      
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" = "NA"))
      temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" = "NA"))

      
    }
  } else {
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_all", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_all", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered", "value" = "NA"))
    
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_no_na_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_no_na_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "corr_of_corr_filtered_no_mut", "value" = "NA"))
    temp.corr_of_corr <- rbind(temp.corr_of_corr, data.frame("method" = "dot_product_filtered_no_mut", "value" = "NA"))
  }
  
  temp.corr_of_corr <- cbind(temp.corr_of_corr, data.frame("mutant" = mutant, "partner" = partner))
  mutant_partner_corr_of_corr <- rbind(mutant_partner_corr_of_corr, temp.corr_of_corr)
}
final.list<-list()
final.list[["cluster"]] <- cluster
final.list[["cluster_number"]] <- cluster_number
final.list[["corr_of_corr_df"]] <- mutant_partner_corr_of_corr

outputfilename <- paste0(Sys.Date(), "_cluster_", cluster_number, "_corr_of_corr.RData")
save(final.list, file = outputfilename)
  