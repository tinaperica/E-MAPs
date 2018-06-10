
options(stringsAsFactors = F)
task_n <- 1
#task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
#inputfilename <- paste0("info_for_correlations/", task_n, "_task_info.RData")
inputfilename <- "clustered_correlations/preprocessed_data_for_correlations/1_task_info.RData"
load(inputfilename)
load("2018-03-09_emap_data_for_corr.RData")
#outputfilename <- paste0("20180310_", task_n, "_correlations.RData")
#output_file_path <- file.path("Output/correlations_w_random/", outputfilename)
output_file_path <- file.path("clustered_correlations/correlation_RData/mutants_only.RData") 
correlations_df <- data.frame()

for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  
  ubermap.query1 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query1, ]
  ubermap.query2 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query2, ]
   
  for ( i in seq_along( ubermap[["clusters"]] ) ) {
    
    cluster <- ubermap[["clusters"]][i]
    temp_library_clusters <- ubermap[["library_clusters"]][ubermap[["library_clusters"]][["cluster"]] == cluster, ]
    temp.ubermap.query1 <- ubermap.query1[ubermap.query1[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    
    if (length(temp.ubermap.query1[,1]) > 10) {   #### I changed all of these conditions to > 10, as that is the final condition below, no point in keeping this as just > 0
      
      temp.ubermap.query2 <- ubermap.query2[ubermap.query2[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
      
      if (length(temp.ubermap.query2[,1]) > 10) {
        merged.emap.data <- merge(temp.ubermap.query1[, c("Gene_uniq", "library.ORF","score","random_score", "random_score_2")], 
                              temp.ubermap.query2[, c("Gene_uniq", "library.ORF", "score", "random_score", "random_score_2")], 
                              by = "library.ORF")
        if (length(merged.emap.data[complete.cases(merged.emap.data[, c("score.x", "score.y")]), ][,1]) > 10) {
          correlation <- round(cor(merged.emap.data$score.x, merged.emap.data$score.y, use = "pairwise.complete.obs"), 3)
        } else {
          correlation <- "NA"
        }
        if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score.x", "random_score.y")]), ][,1]) > 10) {
          random_correlation <- round(cor(merged.emap.data$random_score.x, merged.emap.data$random_score.y, use =  "pairwise.complete.obs"), 3)
        } else {
          random_correlation <- "NA"
        }
        if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score_2.x", "random_score_2.y")]), ][,1]) > 10) {
          random_correlation_2 <- round(cor(merged.emap.data$random_score_2.x, merged.emap.data$random_score_2.y, use =  "pairwise.complete.obs"), 3)
        } else {
          random_correlation_2 <- "NA"
        }
        if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score.x", "random_score_2.y")]), ][,1]) > 10) {
          random_correlation_3 <- round(cor(merged.emap.data$random_score.x, merged.emap.data$random_score_2.y, use = "pairwise.complete.obs"), 3)
        } else {
          random_correlation_3 <- "NA"
        }
        if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score_2.x", "random_score.y")]), ][,1]) > 10) {
          random_correlation_4 <- round(cor(merged.emap.data$random_score_2.x, merged.emap.data$random_score.y, use = "pairwise.complete.obs"), 3)
        } else {
          random_correlation_4 <- "NA"
        }
        correlations_df <- rbind(correlations_df, data.frame( 
              "Gene_uniq1" = query1, "Gene_uniq2" = query2, "cluster" = cluster, 
              correlation, random_correlation, random_correlation_2, random_correlation_3, random_correlation_4))
      }
    }
  }
}

save(correlations_df, file = output_file_path)


