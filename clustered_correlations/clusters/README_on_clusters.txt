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