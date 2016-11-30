options(stringsAsFactors = F)
# combine clusters and calculate correlations of correlations 

first_pair <- as.numeric( Sys.getenv( "SGE_TASK_ID" ))
step_size <- as.numeric( Sys.getenv( "SGE_TASK_ID" )) - 1  ## 551
last_pair <- first_pair + step_size

table_of_pairs <- data.frame()
merged_correlations_df <- data.frame()
for (cluster_number in 1:44) {
  load( paste0( "preprocessed_correlations/", "2016-11-26_", cluster_number, "_corr_lst.RData" ) )
  cluster <- corr_lst[["cluster"]]
  table_of_pairs <- rbind(table_of_pairs, corr_lst[["table_of_pairs"]])
  split_cluster <- unlist(strsplit(cluster, "_"))
  GO_cluster <- split_cluster[1]
  temp_corr <- cbind(corr_lst[["correlations"]], data.frame("GO_cluster" = GO_cluster))
  names(temp_corr)[3] <- "subcluster"
  merged_correlations_df <- rbind(merged_correlations_df, temp_corr)
}
mutants <- as.character(unique(table_of_pairs[["mutant"]]))

table_of_pairs <- table_of_pairs[! duplicated(table_of_pairs),]

mutant_partner_corr_of_corr <- data.frame()
for (p in first_pair:last_pair ) {
  mutant <- as.character( table_of_pairs[["mutant"]][p] )
  partner <- as.character( table_of_pairs[["partner"]][p] ) 
  mut_correlations_df <- merged_correlations_df[
    merged_correlations_df[["gene1"]] == mutant, ]
  temp_mut_correlations_df <- merged_correlations_df[
    merged_correlations_df[["gene2"]] == mutant, ]
  names(temp_mut_correlations_df) <- c("gene2", "gene1", "subcluster", "corr", "filter", "GO_cluster")
  temp_mut_correlations_df <- temp_mut_correlations_df[,c("gene1", "gene2", "subcluster", "corr", "filter", "GO_cluster")]
  mut_correlations_df <- rbind( mut_correlations_df, temp_mut_correlations_df)
  mut_correlations_df <- mut_correlations_df[order(
    mut_correlations_df[["gene2"]], mut_correlations_df[["GO_cluster"]], mut_correlations_df[["subcluster"]]), ]
  partner_correlations_df <- merged_correlations_df[
    merged_correlations_df[["gene1"]] == partner, ]
  temp_partner_correlations_df <- merged_correlations_df[
    merged_correlations_df[["gene2"]] == partner, ]
  names(temp_partner_correlations_df) <- c("gene2", "gene1", "subcluster", "corr", "filter", "GO_cluster")
  partner_correlations_df <- rbind( partner_correlations_df, temp_partner_correlations_df)
  partner_correlations_df <- partner_correlations_df[order(
    partner_correlations_df[["gene2"]], partner_correlations_df[["GO_cluster"]], partner_correlations_df[["subcluster"]]), ]
  if ( (nrow(mut_correlations_df) != 0 ) & ( nrow(partner_correlations_df) != 0) ) { 
    merge_for_corr <- merge(mut_correlations_df, partner_correlations_df, by = c("gene2", "subcluster", "GO_cluster"))
    merge_for_corr <- merge_for_corr[! duplicated(merge_for_corr), ]
    merge_for_corr <- merge_for_corr[complete.cases(merge_for_corr),]
    temp.corr_of_corr <- data.frame()
    if ( nrow(merge_for_corr) > 3 ) {
      names(merge_for_corr) <- expression(other_gene, subcluster, GO_cluster, mutant, mutant.corr, filter.mut, partner, partner.corr, filter.partner)
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
final.list[["cluster"]] <- "merged_clusters"
final.list[["cluster_number"]] <- 45
final.list[["corr_of_corr_df"]] <- mutant_partner_corr_of_corr

outputfilename <- paste0(Sys.Date(), "_cluster_", 45, "_", first_pair, "_corr_of_corr.RData")
save(final.list, file = outputfilename)
