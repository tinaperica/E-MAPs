#### run on guybrush or locally, not on the cluster
#### this script preprocesses the correlation data (from many .RData files) and prepares them for correlation of correlations calculation
### preprocessing outputs an .RData file for each complex
#### correlations output looks like this
###  gene1                  gene2         cluster   corr filter
###### A180T          CTER3XFLAG WT COMPASS complex -0.619     NA

options(stringsAsFactors = F)

corr_lst <- list()  ## final, contains cor per cluster as well as vectors of mut, partner, pairs, clustering method info
correlations <- list() ### contains the correlations per cluster
all_genes <- vector()
dir_name <- "clustered_correlations/ppi_complex_clusters_analysis/corr_RData/"  #### directory that contains the correlation data
files_list <- list.files(path = dir_name)
for (i in seq_along( files_list) ) {
  filename = files_list[i]
  file_path <- file.path(dir_name, filename)
  load(file_path)
  #### I don't know why some genes are NA butmakes sense they need to be removed in the next step
  na_corr_df <- correlations_df[ (is.na(correlations_df$gene1) & is.na(correlations_df$gene2)), ]
  correlations_df <- correlations_df[ ! (is.na(correlations_df$gene1) & is.na(correlations_df$gene2)), ]
  if ( nrow(correlations_df) > 0 ) {
    temp_all_genes <- as.character(unique(c(correlations_df$gene1, correlations_df$gene2)))
    all_genes <- unique(append(all_genes, temp_all_genes))
    clean_names <- data.frame(lapply(correlations_df[,1:2], gsub, pattern = "GSP1 - ", replacement = "", perl = T))
    temp_correlations <- cbind(clean_names, data.frame(
      "cluster" = correlations_df$cluster, "corr" = correlations_df$corr, "filter" = correlations_df$filter))
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
#corr_lst[["partners"]] <- gsub(corr_lst[["partners"]], pattern = " - .+$", replacement = "", perl = T)
corr_lst[["table_of_pairs"]] <- expand.grid("mutant" = corr_lst[["mutants"]], "partner" = corr_lst[["partners"]])
clusters <- corr_lst[["clusters"]]
for ( i in seq_along(clusters) ) {
  cluster <- corr_lst[["clusters"]][i]
  corr_lst[["correlations"]] <- correlations[[cluster]][! duplicated(correlations[[cluster]]), ]
  corr_lst[["cluster"]] <- cluster
  save(corr_lst, file = paste0("clustered_correlations/ppi_complex_clusters_analysis/corr_RData/preprocess_corr/", Sys.Date(), "_", i, "_corr_lst.RData"))
}

print( nrow(corr_lst[["table_of_pairs"]]) )
