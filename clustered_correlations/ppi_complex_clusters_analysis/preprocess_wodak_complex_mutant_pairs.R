### this script does all the clusters but for a small subset (203) of pairs per task (times 44 clusters)


options(stringsAsFactors = F)

ubermap <- list()


#### load an orf to gene name index
orf_gene_name_index <- read.delim("orf_gene_GO_sgd_annotation.txt", head = F)
ubermap[["orf_index"]] <- unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)

#### library genes clustered
ubermap[["cluster_method"]] <- "complexes" #
cluster_file <- paste0("clustered_correlations/clusters/Wodak_yeast_complexes/wodak_complex.txt")
ubermap[["library_clusters"]] <- read.delim(cluster_file, head = T)
ubermap[["clusters"]] <- as.character( unique( ubermap[["library_clusters"]][["cluster"]] ))
##############


## load preprocessed ubermap data (preprocessed so that ubermap$Gene is either a Gsp1 mutant or an ORF)
ubermap[["ubermap"]] <- read.delim("preprocessed_ubermap_all_significant.txt", head = T)


all_genes_and_mutants_df <- read.delim("genes_and_mutants_to_test.txt", head = T)
ubermap[["all_genes_and_mutants"]] <- all_genes_and_mutants_df$unique.genes_and_mutants_to_test.
rm(all_genes_and_mutants_df)
ubermap[["pairs"]] <- combn(ubermap[["all_genes_and_mutants"]], 2)


save(ubermap, file = "correlations_pair_preload.RData")
