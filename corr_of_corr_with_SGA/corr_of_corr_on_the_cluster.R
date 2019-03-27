remove(list = ls())
options(stringsAsFactors = F)

task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
taskfilename <- paste0("all_correlations_task_info/", task_n, "_task_info.RData")
outputfilename <- paste0("20190327_", task_n, "_corr_of_corr.RData")
load(taskfilename)
load("pearson_for_corr_of_corr.RData")

output_file_path <- file.path("Output/corr_of_pearson_corr", outputfilename)

corr_of_corr <- data.frame("query_uniq1" = task.info[["pairs"]][1, ],
                              "query_uniq2" = task.info[["pairs"]][2, ], 
                              "corr_of_corr" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double"))

for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  corr_query1 <- corr[corr[["query_allele_name_1"]] == query1, ]
  corr_query2 <- corr[corr[["query_allele_name_1"]] == query2, ]
  merged.data <- merge(corr_query1, corr_query2, by = "query_allele_name_2")
  merged.data <- merged.data[complete.cases(merged.data),]
  pearson_corr_of_corr <- round(cor(merged.data$corr.x, merged_data$corr.y), 4)
  corr_of_corr[p, 3] <- pearson_corr_of_corr
}

save(corr_of_corr, file = output_file_path)

