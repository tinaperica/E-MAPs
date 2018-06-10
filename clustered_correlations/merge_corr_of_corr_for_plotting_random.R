options(stringsAsFactors = F)

corr_of_corr_df <- data.frame()
clusters <- list()
dir_name <- "corr_of_corr_random"
files_list <- list.files(dir_name)
for ( i in seq_along(files_list) ) {
  filename <- files_list[i]
  if ( grepl("_corr_of_corr.RData", filename) ) {
    file_path <- file.path(dir_name, filename)
    load(file_path)
    cluster <- final.list[["cluster"]]
    clusters[[final.list[["cluster_number"]]]] <- cluster
    corr_of_corr_df <- rbind(corr_of_corr_df, cbind(
      final.list[["corr_of_corr_df"]], "cluster_number" = final.list[["cluster_number"]], "cluster" = cluster) )
  }
}
result <- list()
result[["corr_of_corr"]] <- corr_of_corr_df
result[["clusters"]] <- clusters
#check <- final.list[["corr_of_corr_df"]][order(final.list[["corr_of_corr_df"]]$mutant, 
#                           final.list[["corr_of_corr_df"]]$partner, final.list[["corr_of_corr_df"]]$partner),]
#check <- head(check, n=1000)
rm(final.list)
save(result, file = "random_corr_of_corr.RData")
write.table(corr_of_corr_df, file = "random_cluster_all_corr_of_corr.txt", quote = F, sep = "\t", row.names = F)
