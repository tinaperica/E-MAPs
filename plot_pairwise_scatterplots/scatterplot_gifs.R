library(tidyverse)
library(magick)
library(ggrepel)
abs.pmax <- function(x, y) {
  df <- data.frame(x, y)
  df[is.na(df)] <- 0
  pmax <- vector()
  for (i in 1:length(x)) {
    pair <- df[i, c("x", "y")]
    max <- pair[ which.max(abs(pair)) ][1,1]
    pmax <- append(pmax, max)
  }
  return(pmax)
}
load("basic_E-MAP_data/2018-03-09_emap_data_for_corr.RData")
cluster_annotation <- read_tsv("clustered_correlations/clusters/GO_slims_2018-03-05_pearson_complete_clusters.txt")
representative_cluster_list <- read_tsv("choose_clusters/representative_clusters.txt")
representative_clusters <- cluster_annotation %>% 
  filter(cluster %in% representative_cluster_list$representative_clusters)
all <- tibble("cluster" = "whole library", "ORF" = unique(ubermap$ubermap$library.ORF), 
              "gene_name" = unique(ubermap$ubermap$library.gene_name))
representative_clusters <- bind_rows(representative_clusters, all)
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  select(query, library.ORF, score) %>% 
  inner_join(., representative_clusters,  by = c("library.ORF" = "ORF"))

queries_to_plot <- c("GSP1 Y157A", "GSP1 R108L", "GSP1 T34A", "GSP1 T34E", "GSP1 D79S", 
                     "GSP1 R78K", "GSP1 T34G", "GSP1 R108S",
                     "MOG1 YJR074W", "SRM1 YGL097W_TSQ958", "RNA1 YMR235C_TSQ172", "YRB1 YDR002W_TSQ582")
pairs_to_plot <- combn(queries_to_plot, 2) 
clusters <- unique(representative_clusters$cluster)

### make gifs
for (p in seq_along(pairs_to_plot[1,]) ) {
  img <- image_graph(res = 160)
  query1 <- pairs_to_plot[1,p]
  query2 <- pairs_to_plot[2,p]
  df1 <- data %>% filter(query == query1 & cluster == "whole library")
  df2 <- data %>% filter(query == query2 & cluster == "whole library")
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name"))
  df <- df[complete.cases(df),]
  total_pearson <- round(cor(df$score.x, df$score.y, use = "pairwise.complete.obs"), 2)
  for (cl in seq_along(clusters)) {
    if (clusters[cl] != "whole library") {
      df_sub1 <- data %>% filter(query == query1 & cluster == clusters[cl])
      df_sub2 <- data %>% filter(query == query2 & cluster == clusters[cl])
      df_sub <- inner_join(df_sub1, df_sub2, by = c("library.ORF", "cluster", "gene_name"))
      df_sub <- df_sub[complete.cases(df_sub), ]
      sub_pearson <- round(cor(df_sub$score.x, df_sub$score.y, use = "pairwise.complete.obs"), 2)
      df_strong <- df_sub %>% 
        filter(score.x < -5 | score.x > 2.5 | score.y < -5 | score.y > 2.5)
      pl <- ggplot(data = df, mapping = aes(x = score.x, y = score.y)) + 
        geom_point(color = "grey", alpha = 0.7) +
        xlim(min(c(df$score.x, df$score.y), na.rm = T), max(c(df$score.x, df$score.y), na.rm = T)) +
        ylim(min(c(df$score.x, df$score.y), na.rm = T), max(c(df$score.x, df$score.y), na.rm = T)) +
        theme(legend.key.size =  unit(0.25, "in"),
              legend.text = element_text(size = 9)) +
        geom_point(data = df_sub, mapping = aes(score.x, score.y), color = "black", size = 3) +
        geom_point(data = df_strong, mapping = aes(score.x, score.y, color = gene_name), size = 3) +
        ggtitle(clusters[cl], 
                subtitle =  str_c("Pearson corr for whole library = ", total_pearson, "Pearson corr for the GO category = ", sub_pearson, sep = " ")) + 
        xlab(query1) + ylab(query2) + theme_classic()
      print(pl)
    }
  }
  dev.off()
  animation <- image_animate(img, fps = 0.5)
  image_write(animation, str_c("plot_pairwise_scatterplots/", query1, "_", query2,  "_animation.gif", sep = " "))
}


### make facet_wrap plots
x_lim <- c(min(c(data$score), na.rm = T), max(c(data$score), na.rm = T))
plots <- list()
for (p in seq_along(pairs_to_plot[1,]) ) {
  query1 <- pairs_to_plot[1,p]
  query2 <- pairs_to_plot[2,p]
  df1 <- data %>% filter(query == query1)
  df2 <- data %>% filter(query == query2)
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name"))
  correlations <- df %>% 
    group_by(cluster) %>% 
    summarize("correlation" = round(cor(score.x, score.y, use = "pairwise.complete.obs"), 2)) %>% 
    ungroup()
  df <- df %>% inner_join(., correlations, by = "cluster") %>% 
    mutate("title" = str_c(cluster, correlation, sep = " - Pearson corr = ")) %>% 
    mutate("abs_max" = abs.pmax(score.x, score.y), 
        "score_ratio" = ifelse(abs(score.x) < abs(score.y), abs(score.x/score.y), abs(score.y/score.x)))
  df_bg <- df %>% select(-title)
  df_strong <- df %>% filter(abs(abs_max) > 6 & score_ratio > 0.3)
  plots[[p]] <- ggplot(data = df, mapping = aes(x = score.x, y = score.y), color = "black") + 
    geom_point(data = df_bg, aes(score.x, score.y), color = "grey", alpha = 0.7) +
    geom_point() +
    xlim(x_lim[1], x_lim[2]) +
    ylim(x_lim[1], x_lim[2]) +
    xlab(query1) + ylab(query2) + facet_wrap(~ title) + theme_bw() +
    geom_point(data = df_strong, aes(score.x, score.y, color = gene_name))
}
pdf("plot_pairwise_scatterplots/pairwise_correlations_by_GO_clusters.pdf", width = 14)
print(plots)
dev.off()
