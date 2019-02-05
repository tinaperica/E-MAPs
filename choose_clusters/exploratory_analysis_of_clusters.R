library(tidyverse)
library(mclust)
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
asym_pnorm_cropped <- function(data) {  ### return 1 - 2*probability that the value is same or less (for neg scores) or same or greater for pos scores
  prob_weights <- vector()
  for (x in data) {
    if (x < -2) {
      if ( x < 0 ) {
        p <- pnorm(x, mean = 0, sd = 5)
      } else {
        p <- 0.5
      }
    } else {
      if (x > 1) {
        p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      } else {
        p <- 0.5
      }
    }
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}

asym_pnorm <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
    } 
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}
asym_pnorm_signed <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
      prob_weights <- append(prob_weights, -1 * (1 - 2*p))
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      prob_weights <- append(prob_weights, 1 - 2*p)
    } 
  }
  return(prob_weights)
}

w_pearson <- function( x, y) {
  ## pick pairwise observable
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y) # max absolute value E-MAP score of the two
  w_cropped <- asym_pnorm_cropped(max)  ### 0 for E-MAP scores around 0, used for filtering
  w <- asym_pnorm(max)   ## get the weight for the pair
  df <- cbind(df, data.frame(w, w_cropped))
  df_clean <- df[df$w_cropped != 0, ]  ## df_clean is only for filtering
  if (length(df_clean$x) > 3) {
    # x and y weighted means
    mean_x <- sum(x * w) / sum(w)
    mean_y <- sum(y * w) / sum(w)
    # Compute the weighted variance
    vx <- sum( w * (x - mean_x)^2 ) / sum(w)
    vy <- sum( w * (y - mean_y)^2 ) / sum(w)
    # Compute the covariance
    vxy <- sum( w * (x - mean_x) * (y - mean_y) ) / sum(w)
    # Compute the correlation
    w_correlation <- (vxy / sqrt(vx * vy))
  } else {
    w_correlation <- NaN
  }
  return(w_correlation)
}
SoftCosSim <- function(x, y) {
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  df$wx <- asym_pnorm_signed(df$x)
  df$wy <- asym_pnorm_signed(df$y)
  max <- abs.pmax(df$x, df$y)
  w_cropped <- asym_pnorm_cropped(max)
  df <- cbind(df, data.frame(w_cropped))
  df_clean <- df[df$w_cropped != 0, ]
  df$w_sim <- (2 - abs(df$wx - df$wy)) / 2
  x <- df$x
  y <- df$y
  w_sim <- df$w_sim
  if (length(df_clean$x) > 3) {
    soft_cos_sim <-   sum(w_sim * x * y) / ( sqrt( sum(w_sim * x^2) ) * sqrt( sum(w_sim * y^2) ) )
  } else {
    soft_cos_sim <- NaN
  }
  return(soft_cos_sim)
}
compare_wrap <- function(x, y) {
  values <- vector()
  values <- append(values, w_pearson(x, y))
  values <- append(values, SoftCosSim(x, y))
  return(values)
}
load("basic_E-MAP_data/20180630_emap_data_for_corr_all_clusters.RData")

##### correlations between mutants
#load("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/clustered_correlations/correlation_RData/mutants_only.RData")
# mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
#              "H141R","K101R","T34E","R108Y","NTER3XFLAG,WT","CTER3XFLAG,WT","R108G","R108Q",
#              "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
#              "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
#              "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
#              "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
#              "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
#              "N102M","R108D","K129T")
# cluster_annotation <- read_tsv("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/clustered_correlations/clusters/GO_slims_2018-03-05_pearson_complete_clusters.txt", col_names = T)
# mut_correlations_df <- as_tibble(correlations_df) %>%
#   filter(Gene_uniq1 %in% mutants & Gene_uniq2 %in% mutants) %>%
#   mutate(correlation = as.numeric(correlation), random_correlation = as.numeric(random_correlation),
#          random_correlation_2 = as.numeric(random_correlation_2), random_correlation_3 = as.numeric(random_correlation_3),
#          random_correlation_4 = as.numeric(random_correlation_4), "pair" = str_c(Gene_uniq1, " ", Gene_uniq2))
cluster_annotation <- read_tsv("clustered_correlations/clusters/GO_slims_2018-06-30_pearson_complete_clusters.txt")
dir <- "clustered_correlations/mutant_correlations_all_clusters/"
files <- list.files(path = dir)
data <- data.frame()
for (file in files) {
  inpath = file.path(dir, file)
  load(inpath)
  data <- rbind(data, correlations_df)
  rm(correlations_df)
}
data <- as_tibble(data) %>% 
  mutate("pair" = str_c(Gene_uniq1, Gene_uniq2, sep = " "))
clusters_to_keep <- cluster_annotation %>% pull(cluster) %>% unique()
data <- data %>% 
  filter(cluster %in% clusters_to_keep)

mut_correlations <- data %>% 
  select(Gene_uniq1, Gene_uniq2, pair, cluster, "corr" = w_correlation, "random1" = random_w_correlation, "random2" = random_w_correlation_2)
mut_correlations_high <- data %>% 
  select(Gene_uniq1, Gene_uniq2, pair, cluster, "corr" = w_correlation, "random1" = random_high_w_correlation, "random2" = random_high_w_correlation_2)
mut_cos_sim <- as_tibble(data) %>%
  select(Gene_uniq1, Gene_uniq2, pair, cluster, "corr" = soft_cos_sim, "random1" = random_soft_cos_sim, "random2" = random_soft_cos_sim_2)
mut_cos_sim_high <- as_tibble(data) %>%
  select(Gene_uniq1, Gene_uniq2, pair, cluster, "corr" = soft_cos_sim, "random1" = random_high_soft_cos_sim, "random2" = random_high_soft_cos_sim_2)

pairs <- mut_correlations %>% pull(pair) %>% unique  
# gathered_corr_data <- mut_correlations %>%  
#   #select(Gene_uniq1, Gene_uniq2, cluster, w_correlation, soft_cos_sim, random_high_w_correlation, random_high_w_correlation_2, 
#    #      random_high_soft_cos_sim, random_high_soft_cos_sim_2) %>%
#   gather(key = "correlation_type", value = "corr", -cluster, -Gene_uniq1, -Gene_uniq2, -pair)

# corr_diff <- mut_correlations %>% 
#   #mutate("corr_diff1" = correlation - random_correlation, "corr_diff2" = correlation - random_correlation_2, 
#    #      "corr_diff3" = correlation - random_correlation_3, "corr_diff4" = correlation - random_correlation_4) %>%
#   group_by(cluster) %>%
#   summarise(
#     "correlation_whisker_diff" = boxplot(correlation)$stats[5, ] - boxplot(correlation)$stats[1, ], 
#     "random_correlation_whisker_diff" = boxplot(c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)
#     )$stats[5, ] - 
#       boxplot( c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)
#       )$stats[1, ],
#     "correlation_upper_whisker" = boxplot(correlation)$stats[5, ],
#     "random_upper_whisker" = boxplot(c(random_correlation, random_correlation_2, random_correlation_3, 
#                                        random_correlation_4))$stats[5,],
#     "diff_of_whisker_diffs" = correlation_whisker_diff - random_correlation_whisker_diff,
#     "diff_of_upper_whiskers" = correlation_upper_whisker - random_upper_whisker,
#     "correlation_sd" = sd(correlation),
#     "random_sd" = sd(c(random_correlation, random_correlation_2, random_correlation_3, random_correlation_4)),
#     "sd_diff" = correlation_sd - random_sd
#   ) %>%
#   arrange(desc(sd_diff))

distribution_per_cluster <- function(data_tib, name) {
  diff <- data_tib %>% 
    group_by(cluster) %>% 
    summarise(
      "upper_whisker" = boxplot(corr, plot = F)$stats[5, ],
      "random_upper_whisker" = boxplot(c(random1, random2), plot = F)$stats[5,],
      "diff_of_upper_whiskers" = upper_whisker - random_upper_whisker,
      "sd_diff" = sd(corr, na.rm = T) - sd(c(random1, random2), na.rm = T)
    ) %>% 
    ungroup() %>% 
    arrange(desc(sd_diff))
  distribution_differences <- data_tib %>%
    group_by(cluster) %>%
    summarize("pvalue" = t.test(corr, c(random1, random2))$p.value)
  diff <- inner_join(diff, distribution_differences, by = "cluster")
  return(diff)
}
very_good_clusters <- c(
                        "structural constituent of ribosome", "lipid binding", "organelle fission_GO_30_2_mut",
                        "protein phosphorylation", "enzyme regulator activity_GO_30_1_all")
plot_diff <- function(diff, name, soft_whisk_limit, hard_whisk_limit, soft_sd_limit, hard_sd_limit, pvalue_threshold) {
  very_good <- diff %>% 
    filter(cluster %in% very_good_clusters)
  selected <- diff %>% 
    filter(pvalue > pvalue_threshold) %>% 
    filter(
      (sd_diff > soft_sd_limit & diff_of_upper_whiskers > soft_whisk_limit) |
      (sd_diff > hard_sd_limit) |
      (diff_of_upper_whiskers > hard_whisk_limit)
      )
  diff %>%
    ggplot(aes(x = sd_diff, y = diff_of_upper_whiskers), color = "grey") +
    geom_point(color = "grey") + 
    geom_point(data = selected, aes(x = sd_diff, y = diff_of_upper_whiskers), color = "black") +
    geom_point(data = very_good, aes(x = sd_diff, y = diff_of_upper_whiskers, color = cluster)) +
    geom_hline(yintercept = soft_whisk_limit) + geom_vline(xintercept = soft_sd_limit) + 
    geom_hline(yintercept = hard_whisk_limit) + geom_vline(xintercept = hard_sd_limit) + 
    #scale_colour_gradientn(colours = terrain.colors(10)) +
    ggtitle(name)
  #ggsave(str_c("choose_clusters/", name, "_random_to_real_diff_plot.pdf", sep = ""), width = 10)
  return(selected %>% pull(cluster) %>% unique())
}
plot_distributions <- function(data, diff_data, name, selected_clusters) {
  gathered <- data %>% 
    gather(key = "type", value = "corr", -cluster, -Gene_uniq1, -Gene_uniq2, -pair)
  clusters.ordered.by.significance <- diff_data %>% 
    pull(cluster) %>% unique()
  arranged_by_cluster_sig <- gathered %>%
    mutate("category" = factor(cluster, levels = clusters.ordered.by.significance)) %>%
    arrange(category)
  arranged_by_cluster_sig_selected <- arranged_by_cluster_sig %>%
    filter(cluster %in% selected_clusters)
  plot <- arranged_by_cluster_sig_selected %>% 
    ggplot(aes(x = corr, fill = type)) +
    geom_histogram(alpha = 0.5, bins = 60, position = "identity") + facet_wrap(~category) + theme_bw()
  outfile = str_c("choose_clusters/random_vs_real_", name, "_selected_clusters.pdf", sep = "")
  #ggsave(outfile, width = 48, height = 20)
  #ggsave(outfile, width = 16, height = 7)
}
### weighted correlation, all library genes shiffeled for control
#corr_diff <- distribution_per_cluster(mut_correlations, "w_corr_random")
#plot_diff(corr_diff, "w_corr_random", yintercept = 0.25, xintercept = 0.25)
#corr_diff_clusters <- plot_distributions(mut_correlations, corr_diff, "w_corr_random", 
 #       sd_diff_threshold = 0.25, upper_whisker_diff_threshold = 0.25, pvalue_threshold = 0.0)
### weighted correlation, only high scoring library genes shuffled for control
corr_diff_high <- distribution_per_cluster(mut_correlations_high, "w_corr_random_high")
corr_diff_high_clusters <- plot_diff(corr_diff_high, "w_corr_random_high",
                                  soft_whisk_limit = 0.18, hard_whisk_limit = 0.5,
                                 soft_sd_limit = 0.2, hard_sd_limit = 0.23, pvalue_threshold = 0.0)
# corr_diff_high_clusters <- plot_diff(corr_diff_high, "temp", 
#                                   soft_whisk_limit = 0.43, hard_whisk_limit = 0.7,
#                                   soft_sd_limit = 0.3, hard_sd_limit = 0.4, pvalue_threshold = 0.05)
plot_distributions(mut_correlations_high, corr_diff_high,
                          "w_corr_random_high", corr_diff_high_clusters)
#plot_distributions(mut_correlations_high, corr_diff_high,
 #                  "temp", corr_diff_high_clusters)
### soft cosine similarity, all library genes shuffled for control
#sim_diff <- distribution_per_cluster(mut_cos_sim, "soft_cos_sim_random")
#plot_diff(sim_diff, "soft_cos_sim_random", yintercept = 0.4, xintercept = 0.27)
#sim_diff_clusters <- plot_distributions(mut_cos_sim, sim_diff, "soft_cos_sim_random", 
 #         sd_diff_threshold = 0.27, upper_whisker_diff_threshold = 0.4, pvalue_threshold = 0.0)
### soft cosine similarity, only high scoring library genes shuffled for control
sim_diff_high <- distribution_per_cluster(mut_cos_sim_high, "soft_cos_sim_random_high")
sim_diff_high_clusters <- plot_diff(sim_diff_high, "soft_cos_sim_random_high", 
                                    soft_whisk_limit = 0.2, hard_whisk_limit = 0.45,
                                    soft_sd_limit = 0.19, hard_sd_limit = 0.25, pvalue_threshold = 0.0)
plot_distributions(mut_cos_sim_high, sim_diff_high, "soft_cos_sim_random_high", sim_diff_high_clusters)

union_clusters_to_consider <- union(corr_diff_high_clusters, sim_diff_high_clusters)
union_clusters_to_consider <- append(union_clusters_to_consider, c("whole_library", "sig_Gsp1_GI"))
intersect_clusters_to_consider <- intersect(corr_diff_high_clusters, sim_diff_high_clusters)
intersect_clusters_to_consider <- append(intersect_clusters_to_consider, c("whole_library", "sig_Gsp1_GI"))

selected_clusters_to_use <- cluster_annotation %>%
  filter(cluster %in% union_clusters_to_consider)
#write_tsv(selected_clusters_to_use, path = "choose_clusters/2018-07-02_selected_pearson_complete_clusters.txt")

representative_cluster_annotation <- cluster_annotation %>% 
  filter(cluster %in% union_clusters_to_consider)
queries_to_plot <- c("GSP1 Y157A", "GSP1 R108L", "GSP1 T34E", "GSP1 H141R", 
                     "MOG1 YJR074W",
                     "SRM1 YGL097W_TSQ958", "RNA1 YMR235C_TSQ172")
pairs_to_plot <- combn(queries_to_plot, 2) 
representative_clusters <- representative_cluster_annotation %>% 
  pull(cluster) %>%  unique()
representative_clusters <- "sig_Gsp1_GI_15_7_mut"
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  filter(query %in% queries_to_plot) %>% 
  select(query, library.ORF, score) %>% 
  arrange(query, library.ORF) %>% 
  inner_join(., cluster_annotation,  by = c("library.ORF" = "ORF"))
pairwise_df <- tibble("query.x" = character(), "library.ORF" = character(),
                      "score.x" = double(), "cluster" = character(), "gene_name" = character(),
                      "query.y" = character(), "score.y" = double(), "pair" = character()
                      )
for (p in seq_along(pairs_to_plot[1,]) ) {
  query1 <- pairs_to_plot[1,p]
  query2 <- pairs_to_plot[2,p]
  df1 <- data %>% filter(query == query1)
  df2 <- data %>% filter(query == query2)
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name")) %>% 
    mutate("pair" = str_c(query.x, query.y, sep = " "))
  pairwise_df <- bind_rows(pairwise_df, df)
}
plots <- list()
x_lim <- c(min(c(data$score), na.rm = T), max(c(data$score), na.rm = T))

for (i in seq_along(union_clusters_to_consider)) {
  clust <- union_clusters_to_consider[i]
  df_bg <- pairwise_df %>% 
    select(-cluster) %>% unique()
  df <- pairwise_df %>% 
    filter(cluster == clust)
  correlations <- df %>% 
    group_by(pair) %>% 
    summarize("corr" = round(w_pearson(score.x, score.y), 2), 
              "CosSim" = round(SoftCosSim(score.x, score.y), 2)) %>% 
    ungroup()
  df <- df %>% inner_join(., correlations, by = "pair") %>% 
    mutate("title" = str_c(pair, "\ncorr = ", corr, " CosSim  = ", CosSim, sep = "")) %>% 
    mutate("abs_max" = abs.pmax(score.x, score.y), 
           "score_ratio" = ifelse(abs(score.x) < abs(score.y), abs(score.x/score.y), abs(score.y/score.x)))
  df_strong <- df %>% filter((abs(abs_max) > 3 & score_ratio > 0.2) | abs(abs_max) > 4)
  plots[[i]] <- ggplot(data = df, mapping = aes(x = score.x, y = score.y), color = "black") + 
    geom_point() +
    xlim(x_lim[1], x_lim[2]) +
    ylim(x_lim[1], x_lim[2]) +
    facet_wrap(~ title) + theme_bw() +
    xlab("E-MAP score") + ylab("E-MAP score") +
    geom_point(data = df_strong, aes(score.x, score.y, color = gene_name)) +
    ggtitle(clust)
}
pdf("choose_clusters/pairwise_scatterplots_for_selected_cluster_and_representative_queries.pdf", width = 20, height = 12)
print(plots)
dev.off()





representative_cluster_annotation <- cluster_annotation %>% 
  filter(cluster == "transcription_GO_30_4_all")
queries_to_plot <- c("GSP1 Y157A", "GSP1 R108L", "GSP1 R108A", "GSP1 R108S",
                     "GSP1 T34E", "GSP1 T34A", "GSP1 D79S", "GSP1 H141R", "GSP1 H141I", 
                     "GSP1 Y148I",
                     "MOG1 YJR074W",
                     "SRM1 YGL097W_TSQ958", "RNA1 YMR235C_TSQ172")
pairs_to_plot <- combn(queries_to_plot, 2) 
representative_clusters <- representative_cluster_annotation %>% 
  pull(cluster) %>%  unique()
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  filter(query %in% queries_to_plot) %>% 
  select(query, library.ORF, score) %>% 
  arrange(query, library.ORF) %>% 
  inner_join(., cluster_annotation,  by = c("library.ORF" = "ORF"))
pairwise_df <- tibble("query.x" = character(), "library.ORF" = character(),
                      "score.x" = double(), "cluster" = character(), "gene_name" = character(),
                      "query.y" = character(), "score.y" = double(), "pair" = character()
)
for (p in seq_along(pairs_to_plot[1,]) ) {
  query1 <- pairs_to_plot[1,p]
  query2 <- pairs_to_plot[2,p]
  df1 <- data %>% filter(query == query1)
  df2 <- data %>% filter(query == query2)
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name")) %>% 
    mutate("pair" = str_c(query.x, query.y, sep = " "))
  pairwise_df <- bind_rows(pairwise_df, df)
}
plots <- list()
x_lim <- c(min(c(data$score), na.rm = T), max(c(data$score), na.rm = T))

for (i in seq_along(representative_clusters)) {
  clust <- union_clusters_to_consider[i]
  df_bg <- pairwise_df %>% 
    select(-cluster) %>% unique()
  df <- pairwise_df %>% 
    filter(cluster == clust)
  correlations <- df %>% 
    group_by(pair) %>% 
    summarize("corr" = round(w_pearson(score.x, score.y), 2), 
              "CosSim" = round(SoftCosSim(score.x, score.y), 2)) %>% 
    ungroup()
  df <- df %>% inner_join(., correlations, by = "pair") %>% 
    mutate("title" = str_c(pair, "\ncorr = ", corr, " CosSim  = ", CosSim, sep = "")) %>% 
    mutate("abs_max" = abs.pmax(score.x, score.y), 
           "score_ratio" = ifelse(abs(score.x) < abs(score.y), abs(score.x/score.y), abs(score.y/score.x)))
  df_strong <- df %>% filter((abs(abs_max) > 4 & score_ratio > 0.3) | abs(abs_max) > 5)
  plots[[i]] <- ggplot(data = df, mapping = aes(x = score.x, y = score.y), color = "black") + 
    geom_point() +
    xlim(x_lim[1], x_lim[2]) +
    ylim(x_lim[1], x_lim[2]) +
    facet_wrap(~ title) + theme_bw() +
    xlab("E-MAP score") + ylab("E-MAP score") +
    geom_point(data = df_strong, aes(score.x, score.y, color = gene_name)) +
    ggtitle(clust)
}
pdf("choose_clusters/pairwise_scatterplots_for_best_cluster_and_representative_queries.pdf", width = 36, height = 24)
print(plots)
dev.off()









# check which clusters distinguish between Y157A and R108L
single_mut_pair_check <- mut_correlations %>% 
  filter(pair == "R108L Y157A") %>% 
  rowwise() %>% 
  mutate("corr_diff" = mean(c(corr - random1, corr - random2), na.rm = T)) %>% 
  arrange(desc(abs(corr_diff))) %>% 
  slice(1:648) %>% 
  mutate("order_n" = 1:648)
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  filter(query %in% c("GSP1 R108L", "GSP1 Y157A")) %>% 
  select(query, library.ORF, score) %>% 
  inner_join(., cluster_annotation, by = c("library.ORF" = "ORF"))
order_numbers_to_plot <- seq(1, 648, 12)
plots <- list()
for (i in seq_along(order_numbers_to_plot)) {
  n <- order_numbers_to_plot[i]
  ns <- seq(n, n+11, 1)
  clusters_to_plot <- single_mut_pair_check %>% 
    filter(order_n %in% ns) %>% 
    pull(cluster) %>% unique()
  #temp <- cluster_annotation %>% 
   # filter(cluster %in% clusters_to_plot)
  df1 <- data %>% filter(query == "GSP1 Y157A")
  df2 <- data %>% filter(query == "GSP1 R108L")
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name"))
  df_bg <- df %>% 
    select(-cluster) %>% unique()
  df <- df %>% 
    filter(cluster %in% clusters_to_plot)
  correlations <- df %>% 
    group_by(cluster) %>% 
    summarize("corr" = round(w_pearson(score.x, score.y), 2), 
              "CosSim" = round(SoftCosSim(score.x, score.y), 2)) %>% 
    ungroup()
  df <- df %>% inner_join(., correlations, by = "cluster") %>% 
    mutate("title" = str_c(cluster, "\ncorr = ", corr, " CosSim  = ", CosSim, sep = "")) %>% 
    mutate("abs_max" = abs.pmax(score.x, score.y), 
           "score_ratio" = ifelse(abs(score.x) < abs(score.y), abs(score.x/score.y), abs(score.y/score.x)))
  #df_bg <- df %>% select(-title)
  df_strong <- df %>% filter((abs(abs_max) > 3 & score_ratio > 0.2) | abs(abs_max) > 5)
  plots[[i]] <- ggplot(data = df, mapping = aes(x = score.x, y = score.y), color = "black") + 
    geom_point(data = df_bg, aes(score.x, score.y), color = "grey", alpha = 0.5) +
    geom_point() +
    xlim(x_lim[1], x_lim[2]) +
    ylim(x_lim[1], x_lim[2]) +
    xlab("Y157A") + ylab("R108L") + facet_wrap(~ title) + theme_bw() +
    geom_point(data = df_strong, aes(score.x, score.y, color = gene_name))
}
pdf("choose_clusters/pairwise_scatterplots_R108L_and_Y157A.pdf", width = 14, height = 9.3)
print(plots)
dev.off()

### from these scatterplots I can see that these clusters have more gi with Y157A and less with R108L
clusters_for_Y157A <- c("lipid metabolic process_GO_15_3_all", "lipids_GO_30_1_all",
                        "mRNA binding_GO_15_1_all", "RNA binding_GO_15_1_all", "transcription and mRNA processing_GO_15_12_mut",
                        "RNA binding_GO_30_1_all", "mRNA binding_GO_15_1_mut",
                        "protein modification by small protein conjugation or removal_GO_15_1_all", "regulation of transport",
                        "Golgi", "Golgi and ER_GO_15_10_all", "Golgi and ER_GO_30_4_mut", "Golgi and ER_GO_30_4_all",
                        "sig_Gsp1_GI_15_13_mut", "whole_library_15_33_mut")
#clusters_for_Y157A_annotation <- cluster_annotation %>% 
 # filter(cluster %in% clusters_for_Y157A)
data <- as_tibble(ubermap$ubermap) %>% 
  ungroup() %>% 
  mutate("query" = str_c(Gene.gene_name, Gene_uniq, sep = " ")) %>% 
  filter(query %in% queries_to_plot) %>% 
  select(query, library.ORF, score) %>% 
  arrange(query, library.ORF) %>% 
  inner_join(., cluster_annotation,  by = c("library.ORF" = "ORF"))
x_lim <- c(min(c(data$score), na.rm = T), max(c(data$score), na.rm = T))
plots <- list()
for (p in seq_along(pairs_to_plot[1,]) ) {
  query1 <- pairs_to_plot[1,p]
  query2 <- pairs_to_plot[2,p]
  df1 <- data %>% filter(query == query1)
  df2 <- data %>% filter(query == query2)
  df <- inner_join(df1, df2, by = c("library.ORF", "cluster", "gene_name"))
  df_bg <- df %>% 
    select(-cluster) %>% unique()
  df <- df %>% 
    filter(cluster %in% clusters_for_Y157A)
  correlations <- df %>% 
    group_by(cluster) %>% 
    summarize("corr" = round(w_pearson(score.x, score.y), 2), 
              "CosSim" = round(SoftCosSim(score.x, score.y), 2)) %>% 
    ungroup()
  df <- df %>% inner_join(., correlations, by = "cluster") %>% 
    mutate("title" = str_c(cluster, "\ncorr = ", corr, " CosSim  = ", CosSim, sep = "")) %>% 
    mutate("abs_max" = abs.pmax(score.x, score.y), 
           "score_ratio" = ifelse(abs(score.x) < abs(score.y), abs(score.x/score.y), abs(score.y/score.x)))
  df_strong <- df %>% filter((abs(abs_max) > 2 & score_ratio > 0.3) | abs(abs_max) > 5)
  plots[[p]] <- ggplot(data = df, mapping = aes(x = score.x, y = score.y), color = "black") + 
    geom_point(data = df_bg, aes(score.x, score.y), color = "grey", alpha = 0.5) +
    geom_point() +
    xlim(x_lim[1], x_lim[2]) +
    ylim(x_lim[1], x_lim[2]) +
    xlab(query1) + ylab(query2) + facet_wrap(~ title) + theme_bw() +
    geom_point(data = df_strong, aes(score.x, score.y, color = gene_name))
}
pdf("choose_clusters/pairwise_scatterplots_clusters_that_separate_Y157A.pdf", width = 21, height = 16)
print(plots)
dev.off()






#### how many library genes are covered by these clusters?
lib_genes_considered <- cluster_annotation %>% 
  filter(cluster %in% intersect_clusters_to_consider) %>% 
  arrange(gene_name) %>% 
  pull(gene_name) %>% unique()
length(lib_genes_considered)
##### what are the sizes of the clusters
(cluster_count <- cluster_annotation %>%
    filter(cluster %in% intersect_clusters_to_consider) %>% 
    group_by(cluster) %>%
    summarize("count" = n()) %>% 
    arrange(desc(count)))
tail(cluster_count)
ggplot(cluster_count, aes(x = count)) + geom_histogram(binwidth = 5) + ggtitle("Distribution of cluster sizes")
##### number of clusters per library gene
(n_clusters_per_protein <- cluster_annotation %>%
    filter(cluster %in% intersect_clusters_to_consider) %>% 
    group_by(ORF, gene_name) %>%
    summarize("count" = n()) %>%
    arrange(desc(count)))
n_clusters_per_protein %>% ggplot(aes(x = count)) + geom_histogram(bins = 50) + ggtitle("Number of different clusters per gene")

### how much do the selected clusters overlap?
overlap_between_clusters <- tibble("pair_of_clusters" = character(), 
                                   "overlapping_genes" = character(), cl1_count = double(), cl2_count = double())
cluster_pairs <- combn(intersect_clusters_to_consider, 2)
for (i in 1:(length(cluster_pairs)/2)) {
  cl1 <- cluster_pairs[1, i]
  cl2 <- cluster_pairs[2, i]
  overlapping_genes <- intersect(cluster_annotation$gene_name[cluster_annotation$cluster == cl1], 
                cluster_annotation$gene_name[cluster_annotation$cluster == cl2])
  cl1_count <- length(cluster_annotation$gene_name[cluster_annotation$cluster == cl1])
  cl2_count <- length(cluster_annotation$gene_name[cluster_annotation$cluster == cl2])
  overlap_between_clusters <- bind_rows(overlap_between_clusters, 
        tibble("pair_of_clusters" = str_c(cl1, cl2, sep = " "), overlapping_genes, cl1_count, cl2_count))
}
overlap_size <- overlap_between_clusters %>%
  arrange(cl1_count) %>% 
  group_by(pair_of_clusters, cl1_count, cl2_count) %>% 
  summarize("count" = n())
overlap_size %>% ggplot(aes(x = count)) + geom_histogram(binwidth = 20, alpha = 0.5) + 
  ggtitle("gene overlap between selected clusters") +
  geom_histogram(data = cluster_count, aes(x = count), fill = "blue", alpha = 0.5, binwidth = 20)
overlap_size %>%
  filter(count > 15) %>% 
  ggplot(aes(x = count)) + geom_histogram(binwidth = 1, alpha = 0.5) + 
  ggtitle("gene overlap between selected clusters") +
  geom_histogram(data = cluster_count[cluster_count$count < 100,], aes(x = count), fill = "blue", alpha = 0.5, binwidth = 1)






# arranged_by_cluster_sig <- gathered_corr_data %>%
#   mutate("category" = factor(cluster, levels = clusters.ordered.by.significance)) %>%
#   arrange(category)
#plot <- arranged_by_cluster_sig %>%  ggplot(aes(x = corr, color = correlation_type)) + geom_bar() + facet_wrap(~category)
#ggsave("choose_clusters/many_random_vs_real_Gsp1_mut_correlations_by_library_cluster_by_sd.pdf", width = 48, height = 20)

# clusters_to_consider <- corr_diff %>%
#   filter(sd_diff > 0.2 & diff_of_upper_whiskers > 0.4 & pvalue > 0.25) %>%
#   pull(cluster)
# write_tsv(data.frame(clusters_to_consider), path = "~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/clusters_to_consider.txt")
# arranged_by_cluster_sig_selected <- arranged_by_cluster_sig %>%
#   filter(cluster %in% clusters_to_consider)
# plot <- arranged_by_cluster_sig_selected %>% ggplot(aes(x = corr, color = correlation_type), alpha = 0.1) +
#   geom_density() + facet_wrap(~category)
# ggsave("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/many_random_vs_real_Gsp1_mut_correlations_SELECTED_clusters.pdf", width = 48, height = 20)
# 
# plot <- arranged_by_cluster_sig_selected %>% ggplot(aes(x = corr, color = correlation_type), alpha = 0.1) +
#   geom_bar() + facet_wrap(~category)
# ggsave("~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/many_random_vs_real_Gsp1_mut_correlations_SELECTED_clusters_barplots.pdf", width = 48, height = 20)
# 
# 
# selected_cluster_annotation <- cluster_annotation %>%
#   filter(cluster %in% clusters_to_consider)
# write_tsv(selected_cluster_annotation, path = "~/Documents/Gsp1_bioinformatics/E-MAPs/emap_analysis/choose_clusters/selected_cluster_annotations.txt")
# 
# 
# selected_cluster_annotation %>% 
#   group_by(cluster) %>%
#   summarize("count" = n()) %>%
#   ggplot(aes(x = count)) + geom_histogram() + ggtitle("genes per cluster")
# selected_cluster_annotation %>%
#   group_by(ORF) %>%
#   summarise("count" = n()) %>%
#   ggplot(aes(x = count)) + geom_histogram() + ggtitle("clusters per gene")
# genes_in_clusters <- cluster_annotation %>% pull(gene_name) %>% unique()
# length(genes_in_clusters)
# genes_in_selected_clusters <- selected_cluster_annotation %>% pull(gene_name) %>% unique()
# length(genes_in_selected_clusters)
