library(tidyverse)
library(magick)
load("basic_E-MAP_data/2018-03-09_emap_data_for_corr.RData")
cluster_annotation <- read_tsv("clustered_correlations/clusters/GO_slims_2018-03-05_pearson_complete_clusters.txt")
representative_cluster_list <- read_tsv("choose_clusters/representative_clusters.txt")
representative_clusters <- cluster_annotation %>% 
  filter(cluster %in% representative_cluster_list$representative_clusters)
all <- tibble("cluster" = "all", "ORF" = unique(ubermap$ubermap$library.ORF), 
              "gene_name" = unique(ubermap$ubermap$library.gene_name))
representative_clusters <- bind_rows(representative_clusters, all)
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  select(query, library.ORF, score) %>% 
  inner_join(., representative_clusters,  by = c("library.ORF" = "ORF"))

queries_to_plot <- c("GSP1 Y157A", "GSP1 R108L",  "MOG1 YJR074W", 
      "SRM1 YGL097W_TSQ958", "RNA1 YMR235C_TSQ172", "YRB1 YDR002W_TSQ582")
pairs_to_plot <- combn(queries_to_plot, 2) 
clusters <- unique(representative_clusters$cluster)
for (p in seq_along(pairs_to_plot[1,]) ) {
  img <- image_graph(res = 160)
  for (cl in seq_along(clusters)) {
    if (clusters[cl] != "all") {
      query1 <- pairs_to_plot[1,p]
      query2 <- pairs_to_plot[2,p]
      df1 <- data %>% filter(query == query1 & cluster == "all")
      df2 <- data %>% filter(query == query2 & cluster == "all")
      df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name"))
      total_pearson <- round(cor(df$score.x, df$score.y, use = "pairwise.complete.obs"), 2)
      df_sub1 <- data %>% filter(query == query1 & cluster == clusters[cl])
      df_sub2 <- data %>% filter(query == query2 & cluster == clusters[cl])
      df_sub <- inner_join(df_sub1, df_sub2, by = c("library.ORF", "cluster", "gene_name"))
      sub_pearson <- round(cor(df_sub$score.x, df_sub$score.y, use = "pairwise.complete.obs"), 2)
      pl <- ggplot(data = df, mapping = aes(x = score.x, y = score.y)) + 
        geom_point(color = "grey", alpha = 0.7) +
        geom_point(data = df_sub, mapping = aes(score.x, score.y), color = "black", size = 3) +
        ggtitle(clusters[cl], 
        subtitle =  str_c("Pearson all = ", total_pearson, "Pearson GO category = ", sub_pearson, sep = " ")) + 
        xlab(query1) + ylab(query2) + theme_classic()
      print(pl)
    }
  }
  dev.off()
  animation <- image_animate(img, fps = 0.5)
  image_write(animation, paste(query1, "_", query2,  "_animation.gif"))
}




load("clustered_correlations/correlations_pair_preload.RData")
library_clusters <- ubermap[["library_clusters"]]
saveRDS(library_clusters, "plot_pairwise_scatterplots/data/library_clusters.rds")
clusters <- as.character(unique(library_clusters[["cluster"]]))
cluster_GO_df <- data.frame()
for (i in seq_along(clusters) ) {
  cluster <- clusters[i]
  split_clusters <- unlist(strsplit(cluster, "_"))
  cluster_GO_df <- rbind(cluster_GO_df, data.frame(cluster, "GO_cluster" = split_clusters[1]))
}
saveRDS(cluster_GO_df, "plot_pairwise_scatterplots/data/cluster_GO_df.rds")
GO_clusters <- as.character(unique(cluster_GO_df[["GO_cluster"]]))
#GO_clusters <- c("nuclear transport", "ribosome", "transcription and mRNA processing", "chromatin", "cell cycle", "Golgi and ER", "cytoskeleton")                    
GO_cluster_cols <- list(GO_clusters, colors[1:length(GO_clusters)])
ubermap <- ubermap[["ubermap"]]
saveRDS(ubermap, file = "plot_pairwise_scatterplots/ubermap.rds")
Gene_uniq <- as.character(unique(ubermap$Gene_uniq))
saveRDS(Gene_uniq, file = "plot_pairwise_scatterplots/data/unique_gene_names_from_ubermap.rds")
mutants_to_plot <- "T34A T34Q T34E R108L R108Y D79S D79A YJR074W YGL086W YGR081C"
mutants_to_plot <- unlist(strsplit(mutants_to_plot, split = " "))
mutants_to_plot <- append(mutants_to_plot, c("CTER3XFLAG WT", "NTER3XFLAG WT"))
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
  if (mut1 == "YGR081C") { mut1 <- "SLX9"} 
  if (mut2 == "YGR081C") { mut2 <- "SLX9"} 
  correlation <- cor(merge_to_plot[["score1"]], merge_to_plot[["score2"]], use = "complete.obs")
  plot(merge_to_plot[["score1"]], merge_to_plot[["score2"]], 
       xlab = mut1, ylab = mut2, xlim = c(-20,6), ylim = c(-20, 6), pch = 19, col = "gray")
  legend("top", legend = round(correlation,2))
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
