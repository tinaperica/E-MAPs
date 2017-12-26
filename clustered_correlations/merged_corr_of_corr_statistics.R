setwd("clustered_correlations/")

options(stringsAsFactors = F)

corr_of_corr_dir <- "corr_of_corr/"
files_list <- list.files(corr_of_corr_dir)
merged_corr_of_corr <- data.frame()
for (i in seq_along(files_list)) {
  file <- files_list[i]
  if (grepl("2016-11-29_cluster_45", file)) { ### 45 means it's the merged 44 clusters
    file_path <- file.path(corr_of_corr_dir, file)
    load(file_path)
    merged_corr_of_corr <- rbind(merged_corr_of_corr, final.list[["corr_of_corr_df"]])
  }
}
rm(final.list)
total_n_counts <- with(merged_corr_of_corr[merged_corr_of_corr[["method"]] == "corr_of_corr_all",],
                       aggregate(partner, by=list(mutant=mutant), length))
pdf("Number_of_high_corr_of_corr_per_mutant.pdf", width = 10)
par(mar=c(7,4.1,4.1,2.1))
corr_of_corr_methods <- as.character(unique(merged_corr_of_corr[["method"]]))
for (m in seq_along(corr_of_corr_methods)) {
  method_corr_of_corr_subset <- merged_corr_of_corr[
    merged_corr_of_corr[["method"]] == corr_of_corr_methods[m],]
  if ( length(method_corr_of_corr_subset[["value"]]) > 0) { 
    high_correlation_threshold <- quantile(as.numeric(method_corr_of_corr_subset[["value"]]), prob = .9, na.rm = T)
    filtered_corr_of_corr <- method_corr_of_corr_subset[
        as.numeric(method_corr_of_corr_subset[["value"]]) > high_correlation_threshold, ]
    counts <- with(filtered_corr_of_corr, aggregate(partner, by = list(mutant = mutant), length))
    counts <- counts[order(counts$x, decreasing = T), ]
    high_corr_of_corr_barplot <- barplot(counts$x, col = "royalblue3", 
            main = paste0(corr_of_corr_methods[m], " - number of 10th percentile corr-of-corr interactions (Ntotal = 551)"))
    axis(1, labels = as.character(counts$mutant), at = high_corr_of_corr_barplot, las = 2, cex.axis = 0.7)
  }
}
dev.off()




