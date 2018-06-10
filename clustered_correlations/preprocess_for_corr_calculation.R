library(tidyverse)

# #### inputfiles
# #### SGD index
sgd_index <- "orf_gene_GO_sgd_annotation.txt"
# ### genes and mutants to test
to_test <- "genes_and_mutants_to_test_100merge.txt"
# ### ubermap file
# ### preprocessed ubermap data (preprocessed so that ubermap$Gene is either a Gsp1 mutant or an ORF)
ubermap_input <- "preprocessed_ubermap_all_significant.txt"

#### for random simulated clusters
cluster_method <- "random_clusters"
# ## cluster file
cluster_file <- "clustered_correlations/clusters/simulated_cluster_resampling_control/simulated_cluster_3.txt"
# #### output file
outputfile <- "clustered_correlations/20180122_random_control_correlations_pair_preload.RData"

### for real clusters
#cluster_method <- "GO_slims_pearson_complete"
# #### inputfiles
# ## cluster file
#cluster_file <- "clustered_correlations/clusters/2018-01-17_GO_slims_pearson_complete_clusters.txt"
# #### output file
#outputfile <- "clustered_correlations/20180122_GO_slims_pearson_complete_correlations_pair_preload.RData"

lim.points <- c(-3, 2)  ### this is the threshold for removing the library genes that have significant correlation with wt Gsp1-NAT construct
ubermap <- list()


#### load an orf to gene name index
orf_gene_name_index <- read_delim(sgd_index, "\t", col_names = F)
ubermap[["orf_index"]] <- unique(data.frame("orf" = orf_gene_name_index$X1, "gene_name" = orf_gene_name_index$X2))
### use unique for making orf_index, because GO category term is removed so the remaining table is repetitive 
rm(orf_gene_name_index)

#### library genes clustered
ubermap[["cluster_method"]] <- cluster_method #
ubermap[["library_clusters"]] <- read_delim(cluster_file, delim = "\t", col_names = T)
ubermap[["clusters"]] <- as.character( unique( ubermap[["library_clusters"]][["cluster"]] ))
##############

all_genes_and_mutants_df <- read_delim(to_test, "\t", col_names = T)
names(all_genes_and_mutants_df) <- c("to_test")
ubermap[["all_genes_and_mutants"]] <- all_genes_and_mutants_df$to_test
rm(all_genes_and_mutants_df)

ubermap.data <- read_delim(ubermap_input, "\t")
ubermap.wt <- ubermap.data[grepl( "GSP1-NAT", ubermap.data[["Gene_uniq"]] ), ]
ubermap.wt <- ubermap.wt[findInterval(ubermap.wt[["score"]], lim.points) != 1,]

library_genes_to_remove_due_to_high_score_with_wt <- as.character(unique(ubermap.wt[["library"]]))
ubermap.data <- ubermap.data[! ubermap.data[["library"]] %in% library_genes_to_remove_due_to_high_score_with_wt, ]
ubermap[["ubermap"]] <- ubermap.data

ubermap[["pairs"]] <- combn(ubermap[["all_genes_and_mutants"]], 2)


save(ubermap, file = outputfile)
