options(stringsAsFactors = F)
cluster_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
#cluster_n <- 1
load("mut_gene_corr_of_corr_task_info/1_task_info.RData")  # loads task.info list with task.info[["pairs"]]
#corr_input_file <- paste0("~/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Aug2018_analysis/preprocessed_correlations/weighted_corr_", cluster_n, ".RData")
corr_input_file <- paste0("Output/preprocessed_correlations/weighted_corr_", cluster_n, ".RData")

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

mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T", "NTER3XFLAG WT", "CTER3XFLAG WT")
corr_of_corr <- data.frame()


load(corr_input_file)
#### loads corr_lst with two names corr_lst[["cluster"]] and corr_lst[["correlations"]]
cluster <- corr_lst[["cluster"]]
for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  geneA <- as.character(task.info[["pairs"]][1, p])
  geneB <- as.character(task.info[["pairs"]][2, p])
  #### calculate correlation of correlations for query1 and query2 but not considering correlations between the mutants
  ### don't take into account data where Gene_uniq2 is a Gsp1 mutant
  corrs.geneA <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq1"]] == geneA &
                                               (! corr_lst[["correlations"]][["Gene_uniq2"]] %in% mutants), ]
  temp.corrs.geneA <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq2"]] == geneA &
                                                    (! corr_lst[["correlations"]][["Gene_uniq1"]] %in% mutants), ]
  names(temp.corrs.geneA)[1:2] <- c("Gene_uniq2", "Gene_uniq1") 
  corrs.geneA <- rbind(corrs.geneA, temp.corrs.geneA[, c(2,1,3:8)])
  
  corrs.geneB <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq1"]] == geneB &
                                               (! corr_lst[["correlations"]][["Gene_uniq2"]] %in% mutants), ]
  temp.corrs.geneB <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq2"]] == geneB &
                                                    (! corr_lst[["correlations"]][["Gene_uniq1"]] %in% mutants),  ]
  names(temp.corrs.geneB)[1:2] <- c("Gene_uniq2", "Gene_uniq1")
  corrs.geneB <- rbind(corrs.geneB, temp.corrs.geneB[, c(2,1,3:8)])
  merged.corrs <- merge(corrs.geneA, corrs.geneB, by = "Gene_uniq2")
  if (length(merged.corrs[,1]) > 10) {
    cols <- c(3:8, 10:15)
    merged.corrs[, cols] <- apply(merged.corrs[, cols], 2, function(x) as.numeric(x))
    pear <- pearson_and_p_value(merged.corrs$w_correlation.x, merged.corrs$w_correlation.y)
    random_pear <- pearson_and_p_value(merged.corrs$random_w_correlation.x, merged.corrs$random_w_correlation.y)
    high_pear <- pearson_and_p_value(merged.corrs$random_high_w_correlation.x, merged.corrs$random_high_w_correlation.y)
    soft_sim <- pearson_and_p_value(merged.corrs$soft_cos_sim.x, merged.corrs$soft_cos_sim.y)
    random_soft_sim <- pearson_and_p_value(merged.corrs$random_soft_cos_sim.x, merged.corrs$random_soft_cos_sim.y)
    high_soft_sim <- pearson_and_p_value(merged.corrs$random_high_soft_cos_sim.x, merged.corrs$random_high_soft_cos_sim.y)
    values <- data.frame(matrix(data = c(pear, random_pear, high_pear, soft_sim, random_soft_sim, high_soft_sim),
                                byrow = T, nrow = 6, ncol = 3))
    names(values) = c("corr", "p.value", "sample_size")
    corr_of_corr <- rbind(corr_of_corr, data.frame( 
      geneA, geneB, cluster, 
     "corr_of_corr_type" = c("w_pearson_corr", "random_w_pearson_corr", "random_high_w_pearson_corr", 
              "soft_cos_sim_corr", "random_soft_cos_sim_corr", "random_high_soft_cos_sim_corr"),
     values
     ))
  }
}


output_file_path <- paste0("clustered_correlations/corr_of_corr_mut_partner_per_cluster/20180901_correlations_of_correlations_mut_partner_pairs_only_", 
                           cluster_n,  ".RData")
save(corr_of_corr, file = output_file_path)

