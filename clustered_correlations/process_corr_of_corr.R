options(stringsAsFactors = F)
relevant_clusters <- read.delim("clusters_to_consider.txt", head = T)
task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ))
#task_n <- 1
task.info.filename <- paste0("info_for_corr_of_corr_processing/", task_n, "_task_info.RData") 
load(task.info.filename)
pairs <- unique(task.info[, c("geneA", "geneB")])
files <- unique(task.info$filename)
merged_df <- data.frame()
for (i in seq_along(files))  {    ## this loops through files (one per cluster)
  load(file = files[i])
  final_corr_of_corr <- final_corr_of_corr[, 1:6]
  subset.pairs <- merge(pairs, final_corr_of_corr, by = c("geneA", "geneB"))
  merged_df <- rbind(merged_df, subset.pairs)
}
merged_df$pair <- paste(merged_df$geneA, merged_df$geneB, sep = " ")
merged_df$corr <- as.numeric(merged_df$corr)
all_corr_of_corr <- data.frame()
for (p in seq_along(unique(merged_df$pair))) {
  pair <- unique(merged_df$pair)[p]
  merged_pair <- merged_df[merged_df$pair == pair, c(1:6)]
  temp_pvalues <- merged_pair$p.value
  fdr_p.values <- p.adjust(temp_pvalues, method = "fdr")
  all <- data.frame(merged_pair, "FDR" = fdr_p.values)
  all <- all[all$cluster %in% relevant_clusters$clusters_to_consider, c(1:5, 7)]
  corr.to.keep <- all[all$corr_of_corr_type == "corr_of_corr", ] 
  temp_random <- all[all$corr_of_corr_type != "corr_of_corr",]
  avrg_random_corr <- with(temp_random, aggregate(corr, by = list(cluster = cluster), mean))
  names(avrg_random_corr)[2] <- "random_avrg_corr"
  avrg_random_fdr <- with(temp_random, aggregate(FDR, by = list(cluster = cluster), mean))
  avrg_random_corr <- cbind(avrg_random_corr, "random_avrg_FDR" = avrg_random_fdr$x)
  merged <- merge(corr.to.keep, avrg_random_corr, by = "cluster")
  merged$corr_ratio <- round(abs(merged$corr) / abs(merged$random_avrg_corr), 1)
  all_corr_of_corr <- rbind(all_corr_of_corr, merged)
}

all_outputfilename <- paste0("Output/processed_corr_of_corr/", task_n, ".RData")
save(all_corr_of_corr, file = all_outputfilename)
