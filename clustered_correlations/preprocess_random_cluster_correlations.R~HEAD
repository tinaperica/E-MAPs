options(stringsAsFactors = F)

corr_lst <- list()  ## final, contains corr per cluster as well as vectors of mut, partner, pairs, clustering method info
correlations <- list() ### contains the correlations per cluster
all_genes <- vector()
dir_name <- "random_clusters_final_output_RData"  #### run on guybrush
files_list <- list.files(path = dir_name)
for (i in seq_along( files_list) ) {
  filename = files_list[i]
  file_path <- file.path(dir_name, filename)
  load(file_path)
  correlations_df <- correlations_df[ ! (is.na(correlations_df$gene1) & is.na(correlations_df$gene2)), ]
  if ( nrow(correlations_df) > 0 ) {
    temp_all_genes <- as.character(unique(c(correlations_df$gene1, correlations_df$gene2)))
    all_genes <- unique(append(all_genes, temp_all_genes))
    clean_names <- data.frame(lapply(correlations_df[,1:2], gsub, pattern = "GSP1 - ", replacement = "", perl = T))
    temp_correlations <- cbind(clean_names, correlations_df[,4:6])
    temp_correlations$corr <- as.numeric(temp_correlations$corr)
    corr_lst[["clustering_method"]] <- correlations_df[,3][1]
    rm(correlations_df)  
    clusters <- as.character(unique(temp_correlations$cluster))
    for (j in seq_along(clusters)) {
      cluster <- clusters[j]
      correlations[[cluster]] <- rbind( correlations[[cluster]], 
                      temp_correlations[temp_correlations[["cluster"]] == cluster, ] )
    }
    corr_lst[["clusters"]] <- clusters
  }
}

corr_lst[["mutants"]] <- all_genes[ grepl( "GSP1", all_genes )]
corr_lst[["partners"]] <- all_genes[! all_genes %in% corr_lst[["mutants"]] ]
corr_lst[["mutants"]] <- gsub(corr_lst[["mutants"]], pattern = "GSP1 - ", replacement = "", perl = T)
corr_lst[["table_of_pairs"]] <- expand.grid("mutant" = corr_lst[["mutants"]], "partner" = corr_lst[["partners"]])
clusters <- corr_lst[["clusters"]]
for ( i in seq_along(clusters) ) {
  cluster <- corr_lst[["clusters"]][i]
  corr_lst[["correlations"]] <- correlations[[cluster]][! duplicated(correlations[[cluster]]), ]
  corr_lst[["cluster"]] <- cluster
  save(corr_lst, file = paste0(Sys.Date(), "_", i, "_corr_lst.RData"))
}

print( nrow(corr_lst[["table_of_pairs"]]) )
