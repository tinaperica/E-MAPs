library(tidyverse)
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T", "GSP1-NAT", "NTER3XFLAG WT", "CTER3XFLAG WT")
corr_of_corr_path <- "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Sep2018_analysis/correlations/"
(files <- dir(corr_of_corr_path))
files <- files[grepl(files, pattern = ".RData")]
input_counter <- 0
output_counter <- 0
merged_corr <- data.frame()
for (f in seq_along(files)) {
  file_path <- file.path(corr_of_corr_path, files[f])
  load(file_path)
  correlations_df <- as_tibble(correlations_df)
  merged_corr <- rbind(merged_corr, correlations_df)
  correlations_df <- correlations_df %>% 
    select("Gene_uniq1" = Gene_uniq2, "Gene_uniq2" = Gene_uniq1, w_correlation, soft_cos_sim)
  merged_corr <- rbind(merged_corr, correlations_df)
  input_counter <- input_counter + 1
  if (input_counter == 1581) {
    output_counter <- output_counter + 1
    output_file <- str_c("reverse_corr_of_corr/correlations/correlations_", output_counter, ".RData")
    merged_corr <- merged_corr %>% 
      filter(! Gene_uniq2 %in% mutants)
    save(merged_corr, file = output_file)
    rm(merged_corr)
    merged_corr <- data.frame()
    input_counter <- 0
  }
}
output_file <- str_c("reverse_corr_of_corr/correlations/correlations_10.RData")
save(merged_corr, file = output_file)

### now merge all 10 files and redistribute them by cluster
### also add the random correlations
GO_cat_clusters <- read_tsv("reverse_corr_of_corr/clusters/GO_slims_2018-09-24_pearson_complete_query_genes_clusters.txt")
clusters <- GO_cat_clusters %>% 
  filter(! grepl( "all_queries_", cluster)) %>% 
  pull(cluster) %>% unique()
merged_path <- "reverse_corr_of_corr/correlations"
files <- dir(merged_path)
files <- files[grepl(files, pattern = "correlations_.+RData", perl = T)]
for (cl in seq_along(clusters)) {
  clust <- clusters[cl]
  by_cluster_merged_corr <- list()
  by_cluster_merged_corr[["cluster"]] <- clust
  Gene_uniqs <- GO_cat_clusters %>% 
    filter(cluster == clust) %>% 
    pull(Gene_uniq) %>%  unique()
  for (f in seq_along(files)) {
    file_path <- file.path(merged_path, files[f])
    load(file_path)
    set.seed(2013)
    randomized <- merged_corr %>% 
      mutate("w_correlation" = sample(x = w_correlation, size = length(w_correlation))) %>% 
      mutate("soft_cos_sim" = sample(x = soft_cos_sim, size = length(soft_cos_sim)))
    merged_corr <- merged_corr %>% 
      filter(Gene_uniq2 %in% Gene_uniqs)
    randomized <- randomized %>% 
      filter(Gene_uniq2 %in% Gene_uniqs)
    by_cluster_merged_corr[["corr"]] <- rbind(by_cluster_merged_corr[["corr"]], merged_corr)
    by_cluster_merged_corr[["random"]] <- rbind(by_cluster_merged_corr[["random"]], randomized)
    
  }
  set.seed(2017)
  randomized <- by_cluster_merged_corr[["corr"]] %>% 
    mutate("w_correlation" = sample(x = w_correlation, size = length(w_correlation))) %>% 
    mutate("soft_cos_sim" = sample(x = soft_cos_sim, size = length(soft_cos_sim)))
  by_cluster_merged_corr[["random2"]] <- randomized
  outfilename <- str_c("reverse_corr_of_corr/correlations_by_clusters/cluster_", cl, "_correlations.RData")
  save(by_cluster_merged_corr, file = outfilename)
}

