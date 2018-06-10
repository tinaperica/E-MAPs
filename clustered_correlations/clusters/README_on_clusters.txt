Both sets of clusters (by date) were made by the hclust_EMAP_library_genes.R script
(which performs a range of distance and hierarchical clustering methods)

20161021 clusters are made from the prepocessed_ubermap_ubergenes_only.txt files
(the preprocessed files are all made by the preprocess_ubermap_merged_data.R)


20161101 clusters are made by the same script but from the preprocessed_ubermap_ubergenes_only_significant.txt
--> that file should only contain the library genes that have scores outside of the c(-3,2) range for at least one of the Gsp1 mutants
--> this way, 328 library genes were discarded before clustering was done

NOTE added on 20170130
For all the downstream analysis I used GO slims clusters based on library genes, using pearson (distance) complete (hclust) clustering method
That means I used clusters in this file 2016-11-01_GO_slims_pearson_complete_clusters.txt
i.e. 2018-01-17_GO_slims_pearson_complete_clusters.txt

20180117
    1. in the hclust_EMAP_library_genes.R script remove the option that any protein that
    has "nucleus" as GO category also falls into the nuclear transport category:
            GO_slims_for_clusters[["nuclear transport"]] <- unique ( GO_slims$GO_Slim_term[
                grep( paste( c("nuclear transport"), collapse = "|"), GO_slims$GO_Slim_term) 
                        ])
    instead of
            GO_slims_for_clusters[["nuclear transport"]] <- unique ( GO_slims$GO_Slim_term[
                    grep( paste( c("nuclear transport", "nucleus"), collapse = "|"), GO_slims$GO_Slim_term) 
                        ])
    Now there are only 35 GO (sub)clusters (there is only one nuclear transport cluster)
    2. there is a new preprocessed_ubermap_ubergenes_only_significant.txt file,
        based on basic_E-MAP_data/20180108_gene_names_merge_w_Ubermap_100.txt so need to remake the clusters
        
20180216
    In hclust_EMAP_library_genes.R don't use significant only ubermap
        (significant only only had library genes that have E-MAP scores outside of the c(-3,2)) for at least 1 mutant
        --- the reason I changed my mind and decided to use all the available genes is that even they will have around zero scores
        with each of the mutants they might significantly contribute to negative correlations
    Also, use a more up-to-date GO_slim file 20180216_go_slim_mapping.tab.txt
    