library(RColorBrewer)
colors <- c(brewer.pal(12, "Set3"), brewer.pal(11, "Spectral"))
load("clustered_correlations/correlations_pair_preload.RData")
library_clusters <- ubermap[["library_clusters"]]
clusters <- as.character(unique(library_clusters[["cluster"]]))
cluster_GO_df <- data.frame()
for (i in seq_along(clusters) ) {
  cluster <- clusters[i]
  split_clusters <- unlist(strsplit(cluster, "_"))
  cluster_GO_df <- rbind(cluster_GO_df, data.frame(cluster, "GO_cluster" = split_clusters[1]))
}
#GO_clusters <- as.character(unique(cluster_GO_df[["GO_cluster"]]))
GO_clusters <- c("nuclear transport", "ribosome", "transcription and mRNA processing", "chromatin", "cell cycle", "Golgi and ER", "cytoskeleton")                    
GO_cluster_cols <- list(GO_clusters, colors[1:length(GO_clusters)])
ubermap <- ubermap[["ubermap"]]
mutants_to_plot <- "T34A T34Q T34E R108L R108Y D79S D79A YJR074W YGL086W"
mutants_to_plot <- unlist(strsplit(mutants_to_plot, split = " "))
pairs_to_plot <- combn(mutants_to_plot, 2) 
pdf("pairwise_mutation_E-MAP_score_plots.pdf")
for (p in seq_along(pairs_to_plot[1,]) ) {
  mut1 <- pairs_to_plot[1,p]
  mut2 <- pairs_to_plot[2,p]
  mut1_ubermap <- ubermap[ubermap[["Gene_uniq"]] == mut1, ][,1:3]
  mut2_ubermap <- ubermap[ubermap[["Gene_uniq"]] == mut2, ][,1:3]
  merge_muts <- merge(mut1_ubermap, mut2_ubermap, by = "library")
  names(merge_muts) <- c("library", "mut1", "score1", "mut2", "score2")
  merge_to_plot <- merge(merge_muts, library_clusters, by.x = "library", by.y = "ORF") 
  merge_to_plot <- merge(merge_to_plot, cluster_GO_df, by = "cluster")
  merge_to_plot <- cbind(merge_to_plot, "col" = "red")
  if (mut1 == "YJR074W") { mut1 <- "MOG1"} 
  if (mut2 == "YJR074W") { mut2 <- "MOG1"}
  if (mut1 == "YGL086W") { mut1 <- "MAD1"} 
  if (mut2 == "YGL086W") { mut2 <- "MAD1"} 
  plot(merge_to_plot[["score1"]], merge_to_plot[["score2"]], 
       xlab = mut1, ylab = mut2, xlim = c(-20,6), ylim = c(-20, 6), pch = 19, col = "gray")
  for (g in seq_along(GO_clusters)) {
    GO_cluster <- GO_clusters[g]
    temp <- merge_to_plot[merge_to_plot[["GO_cluster"]] == GO_cluster, ]
    points(temp$score1, temp$score2, col = GO_cluster_cols[[2]][g], pch = 19)
  }
}
plot(merge_to_plot[["score1"]], merge_to_plot[["score2"]], 
     xlab = "", ylab = "", xlim = c(-15,5), ylim = c(-15, 5), type = "n", axes =F)
legend("center", legend = GO_cluster_cols[[1]], fill = GO_cluster_cols[[2]], bty = "n")
dev.off()