options(stringsAsFactors = F)
pearson_and_p_value <- function(x, y) {
  df <- data.frame(x,y)
  df <- df[complete.cases(df), ]
  if (length(df[,1]) > 20) {
    correlation <- round(cor(df$x, df$y, use = "pairwise.complete.obs"), 3)
    p.value <- cor.test(df$x, df$y)$p.value
  } else {
    correlation <- "NA"
    p.value <- "NA"
  }
  return(c(correlation, p.value))
}
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","NTER3XFLAG WT","CTER3XFLAG WT","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T")
cluster_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
inputfilename <- paste0("info_for_mut_partner_corr_of_corr/1_task_info.RData")
load(inputfilename)  # loads task.info list with task.info[["pairs"]]
corr_of_corr <- data.frame()

corr_input_file <- paste0("Output/preprocessed_correlations/corr_w_random__", cluster_n, ".RData")
#corr_input_file <- paste0("clustered_correlations/preprocessed_correlations/corr_w_random__", i, ".RData")

load(corr_input_file)
#### loads corr_lst with two names corr_lst[["cluster"]] and corr_lst[["correlations"]]
cluster <- corr_lst[["cluster"]]
for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  geneA <- task.info[["pairs"]][1, p]
  geneB <- task.info[["pairs"]][2, p]
  #### calculate correlation of correlations for query1 and query2 but not considering correlations between the mutants
  ### don't take into account data where Gene_uniq2 is a Gsp1 mutant
  corrs.geneA <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq1"]] == geneA &
                                               (! corr_lst[["correlations"]][["Gene_uniq2"]] %in% mutants), ]
  temp.corrs.geneA <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq2"]] == geneA &
                                                    (! corr_lst[["correlations"]][["Gene_uniq1"]] %in% mutants), ]
  names(temp.corrs.geneA)[1:2] <- c("Gene_uniq2", "Gene_uniq1") 
  corrs.geneA <- rbind(corrs.geneA, temp.corrs.geneA[, c(2,1,3:7)])
  
  corrs.geneB <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq1"]] == geneB &
                                               (! corr_lst[["correlations"]][["Gene_uniq2"]] %in% mutants), ]
  temp.corrs.geneB <- corr_lst[["correlations"]][ corr_lst[["correlations"]][["Gene_uniq2"]] == geneB &
                                                    (! corr_lst[["correlations"]][["Gene_uniq1"]] %in% mutants),  ]
  names(temp.corrs.geneB)[1:2] <- c("Gene_uniq2", "Gene_uniq1")
  corrs.geneB <- rbind(corrs.geneB, temp.corrs.geneB[, c(2,1,3:7)])
  merged.corrs <- merge(corrs.geneA, corrs.geneB, by = "Gene_uniq2")
  cols <- c(3:7, 9:13)
  merged.corrs[, cols] <- apply(merged.corrs[, cols], 2, function(x) as.numeric(x))
  real <- pearson_and_p_value(merged.corrs$correlation.x, merged.corrs$correlation.y)
  random <- pearson_and_p_value(merged.corrs$random_correlation.x, merged.corrs$random_correlation.y)
  random_2 <- pearson_and_p_value(merged.corrs$random_correlation_2.x, merged.corrs$random_correlation_2.y)
  random_3 <- pearson_and_p_value(merged.corrs$random_correlation_3.x, merged.corrs$random_correlation_3.y)
  random_4 <- pearson_and_p_value(merged.corrs$random_correlation_4.x, merged.corrs$random_correlation_4.y)
  corr_of_corr <- rbind(corr_of_corr, data.frame( 
    geneA, geneB, "pair" = paste0(geneA, geneB), cluster, "corr_of_corr_type" = c("corr_of_corr", "random_corr_of_corr", "random_corr_of_corr_2", 
                                                                                  "random_corr_of_corr_3", "random_corr_of_corr_4"),
    "corr" = c(real[1], random[1], random_2[1], random_3[1], random_4[1]),
    "p.value" = c(real[2], random[2], random_2[2], random_3[2], random_4[2]) ))
}


temp_pairs <- as.character(unique(corr_of_corr$pair))
corr_types <- c("corr_of_corr", "random_corr_of_corr", "random_corr_of_corr_2", 
                "random_corr_of_corr_3", "random_corr_of_corr_4")
final_corr_of_corr <- data.frame()

for (p in seq_along(temp_pairs)) {
  pair <- temp_pairs[p]
  for(corr in seq_along(corr_types)) {
    temp <- corr_of_corr[corr_of_corr$pair == pair & corr_of_corr$corr_of_corr_type == corr_types[corr], c(1,2,4:7)]
    
    temp_pvalues <- temp$p.value
    bonferroni_p.values <- p.adjust(temp_pvalues, method = "bonferroni")
    fdr_p.values <- p.adjust(temp_pvalues, method = "fdr")
    final_corr_of_corr <- rbind(final_corr_of_corr, data.frame(temp, bonferroni_p.values, fdr_p.values))
  }
}


outputfilename <- paste0("20180404_correlations_of_correlations_mut_partner_pairs_only_", cluster_n, ".RData")
outdir <- paste0("Output/corr_of_corr_mut/")
output_file_path <- file.path(outdir, outputfilename)
save(final_corr_of_corr, file = output_file_path)


