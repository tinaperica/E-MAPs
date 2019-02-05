library(tidyverse)
#library(getopt)
#' spec = matrix(c(
#'   'inputdir', 'i', 1, "character",
#'   'outputfile', 'o', 1, "character"
#'   #'all_vs_all', 'a', 0, "logical"
#' ), byrow = TRUE, ncol = 4);
#' opt = getopt(spec)
#if ( is.null(opt$all_vs_all) ) { opt$all_vs_all = FALSE }
#input_corr_dir <- opt$input
input_corr_dir <- "~/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Aug2018_analysis/correlations/"
#outfilename <- opt$outputfile
outfilename <- "clustered_correlations/preprocessed_correlations/weighted_corr"
# mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
#              "H141R","K101R","T34E","R108Y","NTER3XFLAG,WT","CTER3XFLAG,WT","R108G","R108Q",
#              "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
#              "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
#              "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
#              "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
#              "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
#              "N102M","R108D","K129T", "GSP1-NAT", "NTER3XFLAG WT", "CTER3XFLAG WT")
correlations <- list() ### contains the correlations per cluster
all_genes <- vector()
files_list <- list.files(path = input_corr_dir)
for (i in seq_along( files_list) ) {
  filename = files_list[i]
  file_path <- file.path(input_corr_dir, filename)
  load(file_path)
  correlations_df <- as_tibble(correlations_df)
  #correlations[["NA"]] <- correlations_df %>% filter(is.na(Gene_uniq1))
  #correlations_df <- correlations_df %>%
   # filter( ! ((is.na(correlations_df$Gene_uniq1) & is.na(correlations_df$Gene_uniq2)) ))

  if ( nrow(correlations_df) > 0 ) {
    temp_all_genes <- as.character(unique(c(correlations_df$Gene_uniq1, correlations_df$Gene_uniq2)))
    all_genes <- unique(append(all_genes, temp_all_genes))
    clusters <- as.character(unique(correlations_df$cluster))
    for (j in seq_along(clusters)) {
      clust <- clusters[j]
      temp <- correlations_df %>%
        filter(cluster == clust) %>%
        select(Gene_uniq1, Gene_uniq2, w_correlation, soft_cos_sim, random_w_correlation, random_soft_cos_sim, random_high_w_correlation, random_high_soft_cos_sim)
      correlations[[clust]] <- rbind( correlations[[clust]], temp)
    }
  }
}
corr_lst <- list()  ## final
clusters <- names(correlations)
#partners <- all_genes[! all_genes %in% mutants ]
#if (opt$all_vs_all == FALSE) {
 # pairs <- data.frame(expand.grid("gene_A" = mutants, "gene_B" = partners))
#} else {
  ### make all unique combinations, based on this test code: data.frame(t(combn(letters[1:6], m = 2)))
 # combn.temp <- data.frame(t(combn(c(mutants, partners), m = 2)))
#  names(combn.temp) <- c("gene_A", "gene_B")
#  pairs <- combn.temp
#}
for ( i in seq_along(clusters) ) {
  cluster <- clusters[i]
  corr_lst[["correlations"]] <- correlations[[cluster]][! duplicated(correlations[[cluster]]), ]
  corr_lst[["cluster"]] <- cluster
  save(corr_lst, file = paste0(outfilename, "_", i, ".RData"))
}


