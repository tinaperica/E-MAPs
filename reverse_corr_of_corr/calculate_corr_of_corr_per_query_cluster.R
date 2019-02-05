options(stringsAsFactors = F)
cluster_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
load("mut_ms_hit_corr_of_corr_task_info/1_task_info.RData")  # loads task.info list with task.info[["pairs"]]
corr_input_file <- paste0("Output/correlations/cluster_", cluster_n, "_correlations.RData")
pearson_and_p_value <- function(x, y) {
  df <- data.frame(x,y)
  df <- df[complete.cases(df), ]
  n_points <- length(df[,1])
  if (n_points > 10) {
    correlation <- round(cor(df$x, df$y, use = "pairwise.complete.obs"), 4)
    p.value <- cor.test(df$x, df$y)$p.value
  } else {
    correlation <- "NA"
    p.value <- "NA"
  }
  return(c(correlation, p.value, n_points))
}


corr_of_corr <- data.frame()

load(corr_input_file)
#### loads a list with four objects 
### by_cluster_merged_corr[["cluster"]], by_cluster_merged_corr[["corr"]], by_cluster_merged_corr[["random"]] and by_cluster_merged_corr[["random2"]]
cluster <- by_cluster_merged_corr[["cluster"]]
for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  geneA <- as.character(task.info[["pairs"]][1, p])
  geneB <- as.character(task.info[["pairs"]][2, p])
  
  corrs.geneA <- by_cluster_merged_corr[["corr"]][ by_cluster_merged_corr[["corr"]][["Gene_uniq1"]]  == geneA, ]
  corrs.geneB <- by_cluster_merged_corr[["corr"]][ by_cluster_merged_corr[["corr"]][["Gene_uniq1"]]  == geneB, ]
  merged.corrs <- merge(corrs.geneA, corrs.geneB, by = "Gene_uniq2")
  corrs.geneA <- by_cluster_merged_corr[["random"]][ by_cluster_merged_corr[["random"]][["Gene_uniq1"]]  == geneA, ]
  corrs.geneB <- by_cluster_merged_corr[["random"]][ by_cluster_merged_corr[["random"]][["Gene_uniq1"]]  == geneB, ]
  merged.corrs.random <- merge(corrs.geneA, corrs.geneB, by = "Gene_uniq2")
  corrs.geneA <- by_cluster_merged_corr[["random2"]][ by_cluster_merged_corr[["random2"]][["Gene_uniq1"]]  == geneA, ]
  corrs.geneB <- by_cluster_merged_corr[["random2"]][ by_cluster_merged_corr[["random2"]][["Gene_uniq1"]]  == geneB, ]
  merged.corrs.random2 <- merge(corrs.geneA, corrs.geneB, by = "Gene_uniq2")
  if (length(merged.corrs[complete.cases(merged.corrs),1]) > 10) {
    cols <- c(3:4, 6:7)
    merged.corrs[, cols] <- apply(merged.corrs[, cols], 2, function(x) as.numeric(x))
    merged.corrs.random[, cols] <- apply(merged.corrs.random[, cols], 2, function(x) as.numeric(x))
    merged.corrs.random2[, cols] <- apply(merged.corrs.random2[, cols], 2, function(x) as.numeric(x))
    
    pear <- pearson_and_p_value(merged.corrs$w_correlation.x, merged.corrs$w_correlation.y)
    pear_random <- pearson_and_p_value(merged.corrs.random$w_correlation.x, merged.corrs.random$w_correlation.y)
    pear_random2 <- pearson_and_p_value(merged.corrs.random2$w_correlation.x, merged.corrs.random2$w_correlation.y)
    soft_sim <- pearson_and_p_value(merged.corrs$soft_cos_sim.x, merged.corrs$soft_cos_sim.y)
    soft_sim_random <- pearson_and_p_value(merged.corrs.random$soft_cos_sim.x, merged.corrs.random$soft_cos_sim.y)
    soft_sim_random2 <- pearson_and_p_value(merged.corrs.random2$soft_cos_sim.x, merged.corrs.random2$soft_cos_sim.y)
    values <- data.frame(matrix(data = c(pear, pear_random, pear_random2, soft_sim, soft_sim_random, soft_sim_random2),
                                byrow = T, nrow = 6, ncol = 3))
    names(values) = c("corr", "p.value", "sample_size")
    corr_of_corr <- rbind(corr_of_corr, data.frame( 
      geneA, geneB, cluster, 
     "corr_of_corr_type" = c("w_pearson_corr", "random_w_pearson_corr", "random_2_w_pearson_corr", 
              "soft_cos_sim_corr", "random_soft_cos_sim_corr", "random_2_soft_cos_sim_corr"),
     values
     ))
  }
}


output_file_path <- paste0("Output/corr_of_corr_mut_ms_partner_per_cluster/20180925_correlations_of_correlations_mut_partner_ms_pairs_only_", 
                           cluster_n,  ".RData")
save(corr_of_corr, file = output_file_path)

