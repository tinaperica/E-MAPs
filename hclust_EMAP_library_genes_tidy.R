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
ubermap <- read_tsv("basic_E-MAP_data/preprocessed_ubermap_ubergenes_only.txt", col_names = T) 
### 20180216 - don't use significant only
mut.ubermap <- read_tsv("basic_E-MAP_data/preprocessed_ubermap_mut_only.txt", col_names = T)

#### confirm that all Gene_uniq are unique
(ubermap_count <- ubermap %>% group_by(Gene_uniq, library.ORF) %>% 
    summarise("lib_count" = n()) %>% 
    filter(lib_count > 2) %>%
    unique() %>% pull(Gene_uniq)
)
#
##########################
# more checking that everything is complete and unique
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

#### make a group of library gens (and later subclusters of that cluster using hierarchical clustering)
### from all the library genes that have a score spread of more than 7.5 for Gsp1 mutants
## those genes are in library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt
## and this file is made by the preprocess_ubermap_merged_data_tidyverse.R script
range_genes <- read_tsv("library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt", col_names = T)


GO_slims_for_clusters
#### number of GO slims categories all together:
length((names(GO_slims_for_clusters)))
#####################################################################

#### cluster the library genes
### by pearson correlation of e-map scores
### also have clusters where there was no subclustering of GO slims categories 


# ### TO ADD: subcluster everything not on the ubermap scores but on the scores with the Gsp1 mutants
# 
# make_matrices <- function (emap, suffix) {
#   for (cat in names(GO_slims_for_clusters)) {
#     GO_slims_cat <- GO_slims_for_clusters[[cat]]
#     GO_slims_ORFs <- GO_slims %>% 
#       filter(GO_Slim_term %in% GO_slims_cat) %>%
#       pull("ORF") %>% unique()
#     clust.emap <- emap %>% 
#       filter(library.ORF %in% GO_slims_ORFs) 
#     unique_library_genes_with_cat <- clust.emap %>% 
#       pull(library.gene_name) %>% unique()
#     if (length(unique_library_genes_with_cat) > 15) {
#       cluster <- clust.emap %>% 
#         select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
#         mutate("cluster" = cat) %>% unique() %>% 
#         select(cluster, ORF, gene_name)
#       GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
#       clust_emap_spread <- clust.emap %>% 
#         select(Gene_uniq, library.ORF, score) %>% 
#         spread(Gene_uniq, score)
#       temp_matrix <- as.matrix(clust_emap_spread[, -1])
#       rownames(temp_matrix) <- clust_emap_spread$library.ORF
#       ### add each of these clusters onto the lib.ubermap.matrix so I can subscluster them with hclust later
#       lib.ubermap.matrix[[str_c(cat, suffix, "GO", sep = "_")]] <- temp_matrix
#     }
#   }
# }
# ### clustering functions work on matrices
# lib.ubermap.matrix <- list()  # all the matrices are stored in a single list
# GO_cat_clusters <- tibble("cluster" = character(), "ORF" = character(), "gene_name" = character())
# make_matrices(emap = ubermap, suffix = "all")
# make_matrices(mut.ubermap, "mut")

lib.ubermap.matrix <- list()  # all the matrices are stored in a single list
GO_cat_clusters <- tibble("cluster" = character(), "ORF" = character(), "gene_name" = character())
for (cat in names(GO_slims_for_clusters)) {
  GO_slims_cat <- GO_slims_for_clusters[[cat]]
  GO_slims_ORFs <- GO_slims %>% 
    filter(GO_Slim_term %in% GO_slims_cat) %>%
    pull("ORF") %>% unique()
  clust.emap_all <- ubermap %>% 
    filter(library.ORF %in% GO_slims_ORFs)
  clust.emap_mut <- mut.ubermap %>% 
    filter(library.ORF %in% GO_slims_ORFs)
  unique_library_genes_with_cat_all <- clust.emap_all %>% 
    pull(library.gene_name) %>% unique()
  unique_library_genes_with_cat_mut <- clust.emap_mut %>% 
    pull(library.gene_name) %>% unique()
  if (length(unique_library_genes_with_cat_all) > 15) {
    cluster <- clust.emap_all %>% 
      select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, ORF, gene_name) %>% 
      arrange(cluster, ORF)
    GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
  } else if (length(unique_library_genes_with_cat_mut > 15)) {
    cluster <- clust.emap_mut %>% 
      select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, ORF, gene_name)
    GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
  }
  if (length(unique_library_genes_with_cat_all) > 15) {
    cluster <- clust.emap_all %>% 
      select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, ORF, gene_name) %>% 
      arrange(cluster, ORF)
    #GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    clust_emap_spread <- clust.emap_all %>% 
      select(Gene_uniq, library.ORF, score) %>% 
      spread(Gene_uniq, score)
    temp_matrix <- as.matrix(clust_emap_spread[, -1])
    rownames(temp_matrix) <- clust_emap_spread$library.ORF
    ### add each of these clusters onto the lib.ubermap.matrix so I can subscluster them with hclust later
    lib.ubermap.matrix[[str_c(cat, "GO", sep = "_")]][["all"]] <- temp_matrix
  }
  
  if (length(unique_library_genes_with_cat_mut) > 15) {
    cluster <- clust.emap_mut %>% 
      select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, ORF, gene_name)
    #GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    clust_emap_spread <- clust.emap_mut %>% 
      select(Gene_uniq, library.ORF, score) %>% 
      spread(Gene_uniq, score)
    temp_matrix <- as.matrix(clust_emap_spread[, -1])
    rownames(temp_matrix) <- clust_emap_spread$library.ORF
    ### add each of these clusters onto the lib.ubermap.matrix so I can subscluster them with hclust later
    lib.ubermap.matrix[[str_c(cat, "GO", sep = "_")]][["mut"]] <- temp_matrix
  }
}
(GO_cat_clusters <- GO_cat_clusters %>% 
  arrange(ORF, cluster))
#### make a group of library gens (and later subclusters of that cluster using hierarchical clustering)
### from all the library genes that have a score spread of more than 7.5 for Gsp1 mutants
## those genes are in library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt
## and this file is made by the preprocess_ubermap_merged_data_tidyverse.R script
range_genes <- read_tsv("library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt", col_names = T)
range_genes_cluster <- tibble("cluster" = "sig_Gsp1_GI", "ORF" = range_genes$library.ORF, 
                              "gene_name" = range_genes$library.gene_name) 
### add them as a separate category sig_Gsp1_GI to GO_cat_clusters
GO_cat_clusters <- bind_rows(GO_cat_clusters, range_genes_cluster) 
tail(GO_cat_clusters)
(GO_cat_clusters <- GO_cat_clusters %>% 
  arrange(ORF, cluster))

### add them also to the matrix for hierarchical clustering
clust.ubermap <- ubermap %>% 
  filter(library.ORF %in% range_genes$library.ORF) 
clust_ubermap_spread <- clust.ubermap %>% 
  select(Gene_uniq, library.ORF, score) %>% 
  spread(Gene_uniq, score)
temp_matrix <- as.matrix(clust_ubermap_spread[, -1])
rownames(temp_matrix) <- clust_ubermap_spread$library.ORF
lib.ubermap.matrix[["sig_Gsp1_GI"]][["all"]] <- temp_matrix

clust.ubermap.mut <- mut.ubermap %>% 
  filter(library.ORF %in% range_genes$library.ORF) 
clust_ubermap_spread.mut <- clust.ubermap.mut %>% 
  select(Gene_uniq, library.ORF, score) %>% 
  spread(Gene_uniq, score)
temp_matrix <- as.matrix(clust_ubermap_spread.mut[, -1])
rownames(temp_matrix) <- clust_ubermap_spread.mut$library.ORF
lib.ubermap.matrix[["sig_Gsp1_GI"]][["mut"]] <- temp_matrix


outfilename <- str_c("clustered_correlations/clusters/broad_GO_cat_clusters_", Sys.Date(), ".txt", sep = "")
write_tsv(GO_cat_clusters, path = outfilename)

head(lib.ubermap.matrix[["ribosome_GO"]][["mut"]])[,1:20]
head(lib.ubermap.matrix[["ribosome_GO"]][["all"]])[,1:20]

### add the whole library as a cluster
whole_library_genes_cluster <- ubermap %>% 
  select("ORF" = library.ORF, "gene_name" = library.gene_name) %>% 
  unique() %>% 
  mutate("cluster" = "whole_library")
### add as a separate category "whole_library" to GO_cat_clusters
GO_cat_clusters <- bind_rows(GO_cat_clusters, whole_library_genes_cluster) 
tail(GO_cat_clusters)
(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(ORF, cluster))

### Whole ubermap as a matrix
#### library genes, e-map score only matrix
ubermap_spread <- ubermap %>% 
  select(Gene_uniq, library.ORF, score) %>% 
  spread(Gene_uniq, score)
temp_matrix <- as.matrix(ubermap_spread[, -1])
rownames(temp_matrix) <- ubermap_spread$library.ORF
lib.ubermap.matrix[["whole_library"]][["all"]] <- temp_matrix

ubermap_spread.mut <- mut.ubermap %>% 
  select(Gene_uniq, library.ORF, score) %>% 
  spread(Gene_uniq, score)
temp_matrix <- as.matrix(ubermap_spread.mut[, -1])
rownames(temp_matrix) <- ubermap_spread.mut$library.ORF
lib.ubermap.matrix[["whole_library"]][["mut"]] <- temp_matrix

distance_method <- "pearson"
hclust_method <- "complete"
### this function makes subclusters by cutting the tree using cutreeDynamic, 
### it takes the minimal size of the cluster as an argument
## and also deepSplit (from documentation: "The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method")
make_subclusters <- function(cat, input_matrix, GO_cat_name, dendro, dissim, min_size, split_n) { 
  cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, 
                              minClusterSize = min_size, method = "hybrid", deepSplit = split_n, 
                              pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
  clustering.result <- tibble("subcluster" = cut_hybrid, "ORF" = rownames(input_matrix), "GO" = GO_cat_name)
  unclustered_genes <- clustering.result %>% filter(subcluster == 0)
  clusters <- clustering.result %>% filter(subcluster != 0)
  clusters <- clusters %>% 
    mutate("cluster" = str_c(GO_cat_name, min_size, subcluster, cat, sep = "_")) %>%
    select(cluster, ORF) %>%
    inner_join(orf_index, by = "ORF") %>%
    arrange(desc(cluster), ORF)
  return(clusters)
}
### do the clustering for each of the matrices
pdf(paste0("clustered_correlations/clusters/", Sys.Date(), "_clusters.pdf"), width = 10)
for (GO_cat in names(lib.ubermap.matrix)) {
  temp.mat_all <- lib.ubermap.matrix[[GO_cat]][["all"]]
  temp.mat_mut <- lib.ubermap.matrix[[GO_cat]][["mut"]]
  dissim_all <- Dist(temp.mat_all, method = distance_method)
  dissim_mut <- Dist(temp.mat_mut, method = distance_method)
  dendro_all <- hclust(dissim_all, method = hclust_method)
  dendro_mut <- hclust(dissim_mut, method = hclust_method)
  plot(dendro_all, cex = 0.2, main = paste(GO_cat, distance_method, hclust_method, sep = "_"))
  clusters_all <- tibble("cluster" = character(), "ORF" = character(), "gene_name" = character())
  clusters_all <- make_subclusters(cat = "all", input_matrix = temp.mat_all, GO_cat_name = GO_cat, 
                    dendro = dendro_all, dissim = dissim_all, min_size = 15, split_n = 4)
  clusters_all <- bind_rows(clusters_all, make_subclusters(cat = "all", input_matrix = temp.mat_all, GO_cat_name = GO_cat,
                    dendro = dendro_all, dissim = dissim_all, min_size = 30, split_n = 2))
  clusters_all <- bind_rows(clusters_all, make_subclusters(cat = "mut", input_matrix = temp.mat_mut, GO_cat_name = GO_cat,
                    dendro = dendro_mut, dissim = dissim_mut, min_size = 15, split_n = 2))
  clusters_all <- bind_rows(clusters_all, make_subclusters(cat = "mut", input_matrix = temp.mat_mut, GO_cat_name = GO_cat,
                    dendro = dendro_mut, dissim = dissim_mut, min_size = 30, split_n = 4))
  ### check if any of the subclusters are identical
  # if there are any, obviously discard the redundant one
  if (length(unique(clusters_all$cluster)) > 1) {
    all_cluster_pairs <- combn(unique(clusters_all$cluster), 2)
    for (p in seq_along(all_cluster_pairs[1,]) ) {
      cl1 <- all_cluster_pairs[1,p]
      cl2 <- all_cluster_pairs[2,p]
      cl1_df <- clusters_all %>% filter(cluster == cl1)
      cl2_df <- clusters_all %>% filter(cluster == cl2)
      if ( nrow(cl1_df) == nrow(cl2_df) ) {
        if (identical( cl1_df[, -1], cl2_df[, -1])) {
            clusters_all <- clusters_all %>% 
              filter(cluster != cl2)
              print(cl2)
        }
      }
    }
  }
  GO_cat_clusters <- bind_rows(GO_cat_clusters, clusters_all)
}
dev.off()

(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(ORF, cluster))

GO_cat_clusters %>% pull(cluster) %>% unique() %>% length()

(cluster_count <- GO_cat_clusters %>% 
    group_by(cluster) %>%
    summarize("count" = n()) %>% 
    arrange(desc(count)))
tail(cluster_count)
ggplot(cluster_count, aes(x = count)) + geom_histogram() + ggtitle("Distribution of cluster sizes")

##### number of clusters per library gene
(n_clusters_per_protein <- GO_cat_clusters %>% 
  group_by(ORF, gene_name) %>%
  summarize("count" = n()) %>%
  arrange(desc(count)))
#### how many of the library genes don't belong to any clusters?
n_sorted_lib_orfs <- n_clusters_per_protein %>% pull(ORF) %>% unique() %>% length()
n_lib_orfs <- ubermap %>% pull(library.ORF) %>% unique() %>% length()
(unsorted.library.ORFs <- ubermap %>%
  filter(! library.ORF %in% n_clusters_per_protein$ORF) %>% 
  select(library.ORF, library.gene_name) %>% 
  unique())

## check how many unsorted when I remove the whole_library clusters!
GO_cat_clusters_whole_lib_removed <- GO_cat_clusters %>% 
  filter(! grepl("whole_library", cluster))
(n_clusters_per_protein <- GO_cat_clusters_whole_lib_removed %>% 
    group_by(ORF) %>%
    summarize("count" = n()) %>%
    arrange(desc(count)))
(unsorted.library.ORFs <- ubermap %>%
  filter(! library.ORF %in% n_clusters_per_protein$ORF) %>% 
  select(library.ORF, library.gene_name) %>% 
  unique())
### if I don't consider the whole_library clusters there are 45 unsorted library genes
### most of them seem to be of either unknown or poorly defined function
#### check their distribution of scores with the mutants
## if any have interesting interactions with mutants put make an "interesting_unknown" category out of them
#first check their range of scores with mutants

unsorted_genes_score_range <- mut.ubermap %>%
  filter(Gene.gene_name == "GSP1" & library.gene_name %in% unsorted.library.ORFs$library.gene_name) %>% 
  group_by(library.gene_name, library.ORF) %>% 
  do(setNames(data.frame(t(range(.$score, na.rm = T))), c("score_min", "score_max"))) %>% 
  mutate("score_range" = score_max - score_min) %>% 
  arrange(desc(score_range))
### based on the score range, I'll be very inclusive and take first 21
unknown_interesting <- tibble("cluster" = "unknown_interesting", 
                              "ORF" = unsorted_genes_score_range$library.ORF[1:21],
                              "gene_name" = unsorted_genes_score_range$library.gene_name[1:21]
                              )
GO_cat_clusters <- bind_rows(GO_cat_clusters, unknown_interesting)

(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(ORF, cluster))

(n_clusters_per_protein <- GO_cat_clusters %>% 
    group_by(ORF) %>%
    summarize("count" = n()) %>%
    arrange(desc(count)))
ggplot(n_clusters_per_protein, aes(x = count)) + geom_histogram() + ggtitle("Number of different clusters per gene")

outfilename <- str_c("clustered_correlations/clusters/GO_slims", Sys.Date(), "pearson_complete_clusters.txt", sep = "_")
write_tsv(GO_cat_clusters, path = outfilename)


