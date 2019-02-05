library(tidyverse)
corr_of_corr_path <- "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Aug2018_analysis/corr_of_corr_mut_partner_per_cluster/"
(files <- dir(corr_of_corr_path))
files <- files[grepl(files, pattern = ".RData")]
gene_orf <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F) %>% 
  select("ORF" = X1, "gene_name" = X2) %>% 
  unique()
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T", "GSP1-NAT", "NTER3XFLAG WT", "CTER3XFLAG WT")
ms_hits <- read_tsv("SAINT_MS_hits.txt", col_names = F) %>% 
  rename("ORF" = X1) %>% 
  unique()
add <- c("YGR218W", "YER110C", "YAR002W", "YGL241W", "YDR192C", "YDR335W", "YJR074W", 
         "YGL097W", "YMR235C", "YKL205W", "YHR200W", "YMR308C", "YDR002W")
ms_hits <- bind_rows(ms_hits, tibble("ORF" = add)) %>% 
  unique()
pearson_merged_data_for_gia <- data.frame()
soft_cos_merged_data_for_gia <- data.frame()

#mutants <- mutants[1:7]
#files <- files[1:3]
for (m in seq_along(mutants)) {
  mutant <- mutants[m]
  print(mutant)
  merged_df <- data.frame()
  for (i in seq_along(files))  {    ## this loops through files (one per cluster)
    filepath <- file.path(corr_of_corr_path, files[i])
    load(filepath)
    print(filepath)
    corr_of_corr <- as_tibble(corr_of_corr) %>% 
      filter(geneA == mutant & corr_of_corr_type == "w_pearson_corr") %>% 
      mutate("geneB" = gsub(geneB, pattern = "^[0-9]{1}_", replacement = "", perl = T)) %>% ### this is to get rid of things like 1_YDL155W before orf matching
      separate(col = geneB, into = c("ORF"), sep = "_", remove = F) %>% 
      filter(ORF %in% ms_hits$ORF) %>% 
      inner_join(., gene_orf, by = "ORF")
    merged_df <- rbind(merged_df, corr_of_corr)
  }
  rm(corr_of_corr)
  merged_df <- merged_df %>% 
    mutate("pair" = str_c(geneA, geneB, sep = " ")) %>% 
    #group_by(corr_of_corr_type) %>% 
    mutate("adj_FDR" = p.adjust(p.value, method = "fdr")) %>% 
    #ungroup %>% 
    #group_by(pair, corr_of_corr_type) %>%
    group_by(corr_of_corr_type) %>% 
    mutate("n_cluster_adj_FDR" = p.adjust(p.value, method = "fdr")) %>% 
    ungroup() %>% 
    select(-pair)
  outfilename <- str_c("clustered_correlations/corr_of_corr/", mutant, "_20180915_ms_hits_corr_of_corr.RData")
  save(merged_df, file = outfilename)
  
  # merged_df %>% 
  #   filter(corr_of_corr_type == "w_pearson_corr") %>% 
  #   ggplot(., aes(x = p.value)) + geom_density()
  
  #pearson_high_random <- merged_df %>% 
   # filter(corr_of_corr_type == "random_high_w_pearson_corr") %>% 
  #  select(-corr_of_corr_type, -sample_size, -p.value, -adj_FDR)
  pearson <- merged_df %>% 
    #filter(corr_of_corr_type == "w_pearson_corr") %>% 
    select(-corr_of_corr_type, -sample_size, -p.value, -adj_FDR) %>% 
    #inner_join(., pearson_high_random, by = c("geneA", "geneB", "cluster")) %>% 
    filter(n_cluster_adj_FDR < 0.05)
  pearson_data_for_gia <- pearson %>%
    mutate("gene_cluster" = str_c(geneB, cluster, sep = "_")) %>% 
    select(geneA, geneB, gene_cluster, gene_name, corr)
  pearson_merged_data_for_gia <- bind_rows(pearson_merged_data_for_gia, pearson_data_for_gia)
  #pearson_outfilename <- str_c("clustered_correlations/corr_of_corr/", mutant, "_pearson_20180915_all_corr_of_corr.RData")
  #save(pearson, file = pearson_outfilename)
  #pearson_high <- pearson %>% 
  #  filter(n_cluster_adj_FDR.y > 0.05)
  #pearson_outfilename <- str_c("clustered_correlations/corr_of_corr/", mutant, "_pearson_random_control_20180915_all_corr_of_corr.RData")
  #save(pearson_high, file = pearson_outfilename)
  # soft_cos_sim_random <- merged_df %>% 
  #   filter(corr_of_corr_type == "random_high_soft_cos_sim_corr") %>% 
  #   select(-corr_of_corr_type, -sample_size, -p.value, -adj_FDR)
  # soft_cos_sim <- merged_df %>% 
  #   filter(corr_of_corr_type == "soft_cos_sim_corr") %>% 
  #   select(-corr_of_corr_type, -sample_size, -p.value, -adj_FDR) %>% 
  #   inner_join(., soft_cos_sim_random, by = c("geneA", "geneB", "cluster")) %>% 
  #   filter(n_cluster_adj_FDR.x < 0.05)
  # softcos_outfilename <- str_c("clustered_correlations/corr_of_corr/", mutant, "_soft_cos_sim_20180915_all_corr_of_corr.RData")
  # save(soft_cos_sim, file = softcos_outfilename)
  # soft_cos_data_for_gia <- soft_cos_sim %>%
  #   mutate("gene_cluster" = str_c(geneB, cluster, sep = "_")) %>% 
  #   select(geneA, geneB, gene_cluster, corr.x)
  # soft_cos_merged_data_for_gia <- bind_rows(soft_cos_merged_data_for_gia, soft_cos_data_for_gia)
  # soft_cos_sim_high <- soft_cos_sim %>% 
  #   filter(n_cluster_adj_FDR.y > 0.05)
  # softcos_outfilename <- str_c("clustered_correlations/corr_of_corr/", mutant, "_soft_cos_sim_random_control_20180915_all_corr_of_corr.RData")
  # save(soft_cos_sim_high, file = softcos_outfilename)
}

pearson_merged_data_for_gia <- as_tibble(pearson_merged_data_for_gia) 
#soft_cos_merged_data_for_gia <- as_tibble(soft_cos_merged_data_for_gia)

pearson_data_for_gia <- pearson_merged_data_for_gia %>% 
  select(geneA, gene_cluster, corr) %>%
  unique() %>% 
  #mutate("corr" = ifelse(is.na(corr.x), "", corr.x)) %>% 
  spread(gene_cluster, corr, fill = "")

pearson_data_for_gia <- pearson_data_for_gia[, which(colMeans(is.na(pearson_data_for_gia)) < 0.2)]
write_tsv(pearson_data_for_gia, path = "pearson_corr_of_corr_for_gia.txt", na = "")
pearson_annotation_for_gia <- pearson_merged_data_for_gia %>% 
  #mutate("geneB" = gsub(geneB, pattern = "^[0-9]{1}_", replacement = "", perl = T)) %>% ### this is to get rid of things like 1_YDL155W before orf matching
  #separate(col = geneB, into = c("ORF"), sep = "_", remove = F) %>% 
  #inner_join(., gene_orf, by = "ORF") %>% 
  select(gene_cluster, gene_name) %>%  unique()
write_tsv(pearson_annotation_for_gia, path = "pearson_corr_of_corr_annotation_for_gia.txt", col_names = F)

# soft_cos_data_for_gia <- soft_cos_merged_data_for_gia %>% 
#   select(geneA, gene_cluster, corr.x) %>%
#   unique() %>% 
#   #select(geneA, gene_cluster, corr) %>% 
#   spread(gene_cluster, corr.x, fill = "")
# soft_cos_data_for_gia <- soft_cos_data_for_gia[, which(colMeans(is.na(soft_cos_data_for_gia)) < 0.2)]
# write_tsv(soft_cos_data_for_gia, path = "soft_cos_corr_of_corr_for_gia.txt", na = "")
# soft_cos_annotation_for_gia <- soft_cos_merged_data_for_gia %>% 
#   mutate("geneB" = gsub(geneB, pattern = "^[0-9]{1}_", replacement = "", perl = T)) %>% ### this is to get rid of things like 1_YDL155W before orf matching
#   separate(col = geneB, into = c("ORF"), sep = "_", remove = F) %>% 
#   inner_join(., gene_orf, by = "ORF") %>% 
#   select(gene_cluster, gene_name) %>%  unique()
# write_tsv(soft_cos_annotation_for_gia, path = "soft_cos_corr_of_corr_annotation_for_gia.txt", col_names = F)

# pairs <- merged_df %>% 
#   pull(pair) %>% unique()
# measures <- merged_df %>% 
#   pull(corr_of_corr_type) %>% unique()
# 
# all_corr_of_corr <- data.frame()
# for (p in seq_along(pairs)) {
#   pair <- pairs[p]
#   data_for_pair <- merged_df %>% 
#     filter(pair == pair) %>% 
#     select(-pair) %>% 
#     mutate("FDR" = p.adjust(p.value, method = "fdr"))
#   
#   corr.to.keep <- all[all$corr_of_corr_type == "corr_of_corr", ] 
#   temp_random <- all[all$corr_of_corr_type != "corr_of_corr",]
#   avrg_random_corr <- with(temp_random, aggregate(corr, by = list(cluster = cluster), mean))
#   names(avrg_random_corr)[2] <- "random_avrg_corr"
#   avrg_random_fdr <- with(temp_random, aggregate(FDR, by = list(cluster = cluster), mean))
#   avrg_random_corr <- cbind(avrg_random_corr, "random_avrg_FDR" = avrg_random_fdr$x)
#   merged <- merge(corr.to.keep, avrg_random_corr, by = "cluster")
#   merged$corr_ratio <- round(abs(merged$corr) / abs(merged$random_avrg_corr), 1)
#   all_corr_of_corr <- rbind(all_corr_of_corr, merged)
# }
# 
# all_outputfilename <- paste0("Output/processed_corr_of_corr/", task_n, ".RData")
# save(all_corr_of_corr, file = all_outputfilename)
