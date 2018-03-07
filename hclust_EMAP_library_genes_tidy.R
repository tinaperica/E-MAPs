### this script makes hierarchical clusters of genes based on ubermap e-map scores
### as well as based on GO slims terms

library(fastcluster)  # overwrites the hclust function and makes it faster
library(dynamicTreeCut)   ### decides where to cut the hclust tree to make clusters
library(amap)
library(tidyverse)

### ORF to gene name annotation from SGD
orf_gene_name_index <- read_delim("orf_gene_GO_sgd_annotation.txt", col_names = F, delim = "\t")
orf_index <- unique(tibble("ORF" = orf_gene_name_index$X1, "gene_name" = orf_gene_name_index$X2))
rm(orf_gene_name_index)
#############################

###### inputting ubermap data
ubermap <- read_delim("basic_E-MAP_data/preprocessed_ubermap_ubergenes_only.txt", col_names = T, delim = "\t") ### 20180216 - don't use significant only

#### confirm that all Gene_uniq are unique
(ubermap_count <- ubermap %>% group_by(Gene_uniq, library.ORF) %>% 
    summarise("lib_count" = n()) %>% 
    filter(lib_count > 2) %>%
    unique() %>% pull(Gene_uniq)
)
#
##########################
ubermap_per_lib_orf_count <- ubermap %>% group_by(library.ORF) %>%
  summarise("gene_count" = n())
unique(ubermap_per_lib_orf_count$gene_count)
ubermap_per_gene_count <- ubermap %>% group_by(Gene_uniq) %>%
  summarise("lib_count" = n())
unique(ubermap_per_gene_count$lib_count)
ubermap_per_gene_count %>% 
  filter(lib_count > mean(ubermap_per_gene_count$lib_count)) %>%
  pull(Gene_uniq)


### first make clusters based on GO slims
GO_slims <- read_delim( "clustered_correlations/20180216_go_slim_mapping.tab.txt", col_names = F, delim = "\t")
GO_slims <- GO_slims %>% select("ORF" = X1, "Gene" = X2, "GO_Slim_term" = X5)
#### check for overly general slims
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
  summarize("count" = n()) %>% arrange(desc(count)))
### remove some terms
##############################################
GO_slim_terms_to_remove <- as.character(expression(
  signaling, biological_process, response, "response to chemical", "protein complex biogenesis", other, "not_yet_annotated",
  "ion binding", "cellular_component", "molecular_function", membrane, cytoplasm, nucleus, "structural molecule activity", "ATPase activity"
))
GO_slims <- GO_slims %>% filter( ! (GO_Slim_term %in% GO_slim_terms_to_remove))
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
  summarize("count" = n()) %>% arrange(desc(count)))
tail(GO_slims_count)
ggplot(GO_slims_count, aes(x = count)) + geom_histogram() + ggtitle("Go slims terms / cluster sizes")
GO_slim_terms <- GO_slims %>% pull(GO_Slim_term) %>% unique()
### make each remaining GO_slim term it's own cluster
#### and then make some composite clusters (e.g. combin ER and Golgi, or combine multiple terms into mRNA processing cluster etc.)
##################################################################
## first add all the terms as clusters
GO_slims_for_clusters <- list()
for (i in seq_along(GO_slim_terms)) {
  GO_slims_for_clusters[[GO_slim_terms[i]]] <- GO_slim_terms[i]
}
### then manually select some terms based on which I will make GO based clusters
#### make sure that composite GO_slims_terms have unique names, that don't already exist as GO_slims_terms
##### this approach gives a small number of relativelly general clusters
GO_slims_for_clusters[["ribosomes and translation"]] <- unique ( GO_slims$GO_Slim_term[ 
  grep( paste( c("rRNA", "ribosom", "translation", "tRNA", "snoRNA"), collapse = "|"), GO_slims$GO_Slim_term)
  ])
GO_slims_for_clusters[["transcription and mRNA processing"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mRNA", "RNA splicing", "transcription", "RNA modification", "mRNA processing",
                 "rRNA processing", "tRNA processing", "snoRNA processing"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["transcription"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("transcription"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["RNA processing"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mRNA processing", "splicing", "snoRNA processing"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["Golgi and ER"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("protein lipidation", "protein maturation", "endocytosis", "regulation of transport", "glycosylation",
                 "vesicle organization","endosom", "Golgi", "endoplasmic"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["Golgi"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("Golgi"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["ER"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("endoplasmic reticulum"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["peroxisomes"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("peroxisom"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["vacuoles"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("vacuol"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["mitochondria"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mitochond", "respirat"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["chromatin"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("chromatin","histone", "telomere", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["cytoskeleton and microtubules"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cytoskelet", "cytokinesis", "chromosome segregation", "microtubule"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["cell cycle"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cell cycle", "recombination", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["budding"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("bud", "fission", "sporulation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["lipids"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("lipid"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["nuclear transport and organization"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("nuclear transport", "nucleus"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
GO_slims_for_clusters[["metabolic"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("metabolic", "metabolite"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
#####################
GO_slims_for_clusters
#### number of GO slims categories all together:
length((names(GO_slims_for_clusters)))
#####################################################################

#### cluster the library genes
### by pearson correlation of e-map scores
### also have clusters where there was no subclustering of GO slims categories 

### clustering functions work on matrices
lib.ubermap.matrix <- list()  # all the matrices are stored in a single list

GO_cat_clusters <- tibble("cluster" = character(), "ORF" = character(), "gene_name" = character())
for (cat in names(GO_slims_for_clusters)) {
  GO_slims_cat <- GO_slims_for_clusters[[cat]]
  GO_slims_ORFs <- GO_slims %>% filter(GO_Slim_term %in% GO_slims_cat) %>% pull("ORF") %>% unique()
  clust.ubermap <- ubermap %>% filter(library.ORF %in% GO_slims_ORFs) 
  unique_library_genes_with_cat <- clust.ubermap %>% pull(library.gene_name) %>% unique()
  if (length(unique_library_genes_with_cat) > 30) {
    cluster <- clust.ubermap %>% select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% select(cluster, ORF, gene_name)
    GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    ### add each of these clusters onto the lib.ubermap.matrix so I can subscluster them with hclust later
    Gene_uniq <- clust.ubermap %>% pull(Gene_uniq) %>% unique()
    library.ORF <- clust.ubermap %>% pull(library.ORF) %>% unique()
    lib.ubermap.matrix[[paste0(cat, "_GO")]] <- matrix( clust.ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(Gene_uniq) )
    rownames( lib.ubermap.matrix[[ paste0( cat, "_GO" ) ]] ) <- library.ORF
    colnames( lib.ubermap.matrix[[ paste0( cat, "_GO" ) ]] ) <- Gene_uniq
  }
}
outfilename <- str_c("clustered_correlations/clusters/broad_GO_cat_clusters", Sys.Date(), ".txt", sep = "")
write_delim(GO_cat_clusters, path = outfilename, delim = "\t")
head(lib.ubermap.matrix[["ribosome_GO"]])[,1:20]


### Whole ubermap as a matrix
#### library genes, e-map score only matrix
Gene_uniq <- ubermap %>% pull(Gene_uniq) %>% unique()
library.ORF <- ubermap %>% pull(library.ORF) %>% unique()
lib.ubermap.matrix[["whole_library"]] <- matrix(ubermap$score, byrow = F, nrow = length(library.ORF), ncol = length(Gene_uniq))
rownames(lib.ubermap.matrix[["whole_library"]]) <- library.ORF
colnames(lib.ubermap.matrix[["whole_library"]]) <- Gene_uniq
#head(lib.ubermap.matrix[["library"]])[,1:20]  ### rows are library gene ORFs and columns are unique Gene identifiers (check preprocess_ubermap_merged_data.R for why uniq)


distance_methods <- c("pearson")
hclust_methods <- c("complete")
### do the clustering for each of the matrices
clustering_summary_df <- data.frame()
pdf(paste0("clustered_correlations/clusters/", Sys.Date(), "_clusters.pdf"), width = 10)
for (i in seq_along(distance_methods)) {
  for (j in seq_along(hclust_methods)) {
    GO_slims_subcluster_unclustered_genes <- vector()
    GO_slims_subcluster_unique_clustered_genes <- vector()
    for (mat in names(lib.ubermap.matrix)) {
      temp.mat <- lib.ubermap.matrix[[mat]]
      dist_method <- distance_methods[i]
      hclust_method <- hclust_methods[j]
      dissim <- Dist(temp.mat, method = dist_method)
      dendro <- hclust(dissim, method = hclust_method)
      plot(dendro, cex = 0.2, main = paste(mat, dist_method, hclust_method, sep = "_"))
      cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, minClusterSize = 30, method = "hybrid", deepSplit = 4, pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
      clustering.result <- tibble("subcluster" = cut_hybrid, "ORF" = rownames(temp.mat), "GO" = mat)
      unclustered_genes <- clustering.result %>% filter(subcluster == 0)
      clusters <- clustering.result %>% filter(subcluster != 0)
      clusters <- clusters %>% 
        mutate("cluster" = str_c(GO, subcluster, sep = "_")) %>%
        select(cluster, ORF) %>%
        inner_join(orf_index, by = "ORF") %>%
        arrange(desc(cluster))
      GO_cat_clusters <- bind_rows(GO_cat_clusters, clusters)
    }
    outfilename <- str_c("clustered_correlations/clusters/GO_slims", Sys.Date(), dist_method, hclust_method, "clusters.txt", sep = "_")
    write_delim(GO_cat_clusters, path = outfilename, delim = "\t")
  }
}
dev.off()



(cluster_count <- GO_cat_clusters %>% 
    group_by(cluster) %>%
    summarize("count" = n()) %>% 
    arrange(desc(count)))
tail(cluster_count)
ggplot(cluster_count, aes(x = count)) + geom_histogram() + ggtitle("Distribution of cluster sizes")

(n_clusters_per_protein <- GO_cat_clusters %>% 
  group_by(ORF) %>%
  summarize("count" = n()) %>%
  arrange(desc(count)))

ggplot(n_clusters_per_protein, aes(x = count)) + geom_histogram() + ggtitle("Number of different clusters per gene")

