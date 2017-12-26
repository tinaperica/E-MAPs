library(gplots)
#### this script outputs heatmaps, one for each Gsp1 partner of interest, where columns are point mutants of Gsp1
### and the rows are different clusters based on GO categories (values in the heatmaps are correlations of correlations)
options(stringsAsFactors = F)

ms_partners_df <- read.delim("Gsp1_MS_hits.txt", head = F)
ms_partners <- ms_partners_df$V1
load("corr_of_corr.RData")
corr_of_corr_df <- result[["corr_of_corr"]]
clusters <- result[["clusters"]]
rm(result)

corr_of_corr_df <- data.frame(lapply(corr_of_corr_df, gsub, pattern = " - [A-z0-9]{7}", replacement = "", perl = T))
corr_of_corr_df <- corr_of_corr_df[ order(corr_of_corr_df$mutant, corr_of_corr_df$partner, corr_of_corr_df$method),]
partners.list<- list()
core_modifiers <- c("RNA1", "SRM1", "YRB1", "PSE1", "MOG1")
partners.list[["Core_modifiers"]] <- core_modifiers
selected_measures <- c( "corr_of_corr_no_na_no_mut")

temp <- corr_of_corr_df[ grepl(pattern = paste(core_modifiers, collapse="|"), corr_of_corr_df[["partner"]]), ]
partners_to_plot <- as.character(unique(temp[["partner"]]))
partners_to_plot <- partners_to_plot[! partners_to_plot == "RNA15_TSQ652"]
colors <- rainbow(length(partners_to_plot))
coloring_list <- list()
for (i in seq_along(partners_to_plot)) {
  coloring_list[partners_to_plot[i]] <- colors[i]
}
print_heatmap <- function(df, rows, columns, file_prefix, file_suffix, file_height, file_width) {
  clusters.to.remove <- unique(df$cluster[df$value == "NA"])
  if (length(clusters.to.remove) > 0) {
    df <- df[! grepl(pattern = paste(clusters.to.remove, collapse="|"), df[["cluster"]]), ]
  }
  df <- df[ order(df[[columns]], df[[rows]]), ]
  heatmap_rows <- as.character( unique( df[[rows]] ))
  heatmap_columns <- as.character(unique( df[[columns]] ))
  if ( nrow(df) > 2) {
    correlation.matrix <- matrix( as.numeric(df$value), byrow = T, nrow = length(heatmap_columns), ncol = length(heatmap_rows) )
    correlation.matrix <- t(correlation.matrix)
    colnames(correlation.matrix) <- heatmap_columns
    rownames(correlation.matrix) <- heatmap_rows
    
    cols <- rep('black', nrow(correlation.matrix))
    for (i in seq_along(names(coloring_list))) {
      cols[ grepl(pattern = names(coloring_list)[i], row.names(correlation.matrix)) ] <- unlist(coloring_list[[i]])
    }
    
    filename <- paste0("clustered_correlations/corr_of_corr_plots/per_partner_heatmaps/", Sys.Date(), "_", file_prefix, "_", measure, "_", file_suffix, ".pdf")
    pdf(file = filename, width = file_width, height = file_height)
    heatmap <- heatmap.2( correlation.matrix, scale = "none", dendrogram="both", 
                          trace="none", density.info="none", col = cm.colors, 
                          key.ylab = rows, key.xlab = columns, 
                          main = file_suffix, cexCol=0.9, cexRow = 0.9, keysize=0.9, margin=c(10,20), colRow = cols)
    
    dev.off()
  }
}
for ( m in seq_along(selected_measures) ) {
  measure = selected_measures[m]
  temp_measure.df <- corr_of_corr_df[
    corr_of_corr_df[["method"]] == measure, ]
  #merged_subclusters_corr_of_corr <- list()
  ### this part plots heatmaps per partner (where rows are different clusters)
  for ( p in seq_along( partners_to_plot ) ) {
    partner <- partners_to_plot[[p]]
    temp_measure_partner.df <- temp_measure.df[ temp_measure.df[["partner"]] == partner,]
    go_subs.df <- temp_measure_partner.df[ grepl(pattern = "GO", temp_measure_partner.df[["cluster"]]), ]
    print_heatmap(go_subs.df, "cluster", "mutant", "GO_subclusters", partner, 14, 14)
    go.df <- temp_measure_partner.df[! grepl(pattern = "GO", temp_measure_partner.df[["cluster"]]), ]
    go.df <- go.df[go.df[["cluster"]] != "merged_clusters", ]
    print_heatmap(go.df, "cluster", "mutant", "GO_clusters", partner, 10, 14)
  }
  ## this part concatanates each of the core modifiers with each of the clusters and uses that as a row in the heatmap
  ##################
  ### Get oly rows with either of the partners to plot
  temp_measure_all_partners.df <- temp_measure.df[ temp_measure.df[["partner"]] %in% partners_to_plot, ]
  # now concatanate the partner and the cluster field
  temp_measure_all_partners.df <- cbind( temp_measure_all_partners.df, 
              "unique" = paste(temp_measure_all_partners.df[["partner"]], temp_measure_all_partners.df[["cluster"]]) )
  go_subs_all_partners.df <- temp_measure_all_partners.df[ grepl(pattern = "GO", temp_measure_all_partners.df[["cluster"]]), ]
  print_heatmap(go_subs_all_partners.df, "unique", "mutant", "GO_subclusters", "all_partners", 35, 14)
  go_all_partners.df <- temp_measure_all_partners.df[ ! grepl(pattern = "GO", temp_measure_all_partners.df[["cluster"]]), ]
  go_all_partners.df <- go_all_partners.df[ go_all_partners.df[["cluster"]] != "merged_clusters", ]
  print_heatmap(go_all_partners.df, "unique", "mutant", "GO_clusters", "all_partners", 14, 20)
}








