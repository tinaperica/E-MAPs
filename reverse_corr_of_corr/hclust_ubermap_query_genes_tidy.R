### this script makes hierarchical clusters of query genes based on GO slims terms as well
# as subclusters of those based on genetic interactions (ubermap e-map scores)

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
#ubermap <- read_tsv("basic_E-MAP_data/preprocessed_ubermap_ubergenes_only.txt", col_names = T) 
### 20180216 - don't use significant only
#mut.ubermap <- read_tsv("basic_E-MAP_data/preprocessed_ubermap_mut_only.txt", col_names = T)
ubermap <- read_tsv("basic_E-MAP_data/20180920_preprocessed_ubermap_500_overlap_ubergenes_only.txt")

Gene_uniq_gene_name_index <- ubermap %>% 
  select(Gene_uniq, Gene.gene_name) %>% 
  unique()
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
GO_slims <- read_delim( "reverse_corr_of_corr/20180921_go_slim_mapping.tab.txt", col_names = F, delim = "\t") %>% 
  select("ORF" = X1, "Gene" = X2, "GO_Slim_term" = X5)
#### check for overly general slims
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
  summarize("count" = n()) %>% arrange(desc(count)))
### remove some terms
##############################################
GO_slim_terms_to_remove <- as.character(expression(
  signaling, biological_process, response, "response to chemical", "protein complex biogenesis", other, 
  "not_yet_annotated", "ion binding", "cellular_component", "molecular_function", 
  membrane, cytoplasm, nucleus, "extracellular region",
  "RNA binding", "DNA binding", "plasma membrane", "protein folding", "enzyme binding", "protein binding, bridging",
  "oxidoreductase activity", "enzyme regulator activity", "protein phosphorylation", "kinase activity",
  "conjugation", "ligase activity", "phosphatase activity", "nuclease activity", "peptidase activity", 
  "lyase activity", "structural molecule activity", "ATPase activity", "hydrolase activity", 
  "transferase activity", "methyltransferase activity", "helicase activity", "isomerase activity", 
  "protein transporter activity", "protein alkylatio", "protein acylation", "protein dephosphorylation", 
  "hydrolase activity, acting on glycosyl bonds"
))
GO_slims <- GO_slims %>% filter( ! (GO_Slim_term %in% GO_slim_terms_to_remove))
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
  summarize("count" = n()) %>% arrange(desc(count)))
tail(GO_slims_count)
ggplot(GO_slims_count, aes(x = count)) + geom_histogram() + ggtitle("GO slims terms / cluster sizes")
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
GO_slims_for_clusters[["mitochondria and respiration"]] <- unique ( GO_slims$GO_Slim_term[
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

#### cluster the query genes
### by pearson correlation of e-map scores
### also have clusters where there was no previous subclustering of GO slims categories 


query_genes.ubermap.matrix <- list()  # all the matrices are stored in a single list
GO_cat_clusters <- tibble("cluster" = character(), "Gene_uniq" = character(), "Gene.gene_name" = character())
for (cat in names(GO_slims_for_clusters)) {
  GO_slims_cat <- GO_slims_for_clusters[[cat]]
  GO_slims_ORFs <- GO_slims %>% 
    filter(GO_Slim_term %in% GO_slims_cat) %>%
    pull("ORF") %>% unique()
  clust.emap_queries <- ubermap %>% 
    filter(ORF %in% GO_slims_ORFs)
  clust_unique_queries <- clust.emap_queries %>% 
    pull(Gene_uniq) %>% unique() %>% sort()
  if (length(clust_unique_queries) >= 50) {
    cluster <- clust.emap_queries %>% 
      select(Gene_uniq, Gene.gene_name) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, Gene_uniq, Gene.gene_name) %>% 
      arrange(cluster, Gene_uniq)
    GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    #GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    clust_emap_spread <- clust.emap_queries %>% 
      select(Gene_uniq, library.ORF, score) %>% 
      spread(library.ORF, score)
    temp_matrix <- as.matrix(clust_emap_spread[, -1])
    rownames(temp_matrix) <- clust_emap_spread$Gene_uniq
    ### add each of these clusters onto the query_genes.ubermap.matrix so I can subscluster them with hclust later
    query_genes.ubermap.matrix[[str_c(cat, "GO", sep = "_")]] <- temp_matrix
  }
}

# check that it looks as expected
head(query_genes.ubermap.matrix[["ribosome_GO"]])[,1:20]

outfilename <- str_c("reverse_corr_of_corr/clusters/broad_GO_cat_clusters_", Sys.Date(), ".txt", sep = "")
write_tsv(GO_cat_clusters, path = outfilename)

### add the total of all the query genes as a "cluster"
all_query_genes_cluster <- ubermap %>% 
  select("Gene_uniq" = Gene_uniq, Gene.gene_name) %>% 
  unique() %>% 
  arrange(Gene.gene_name) %>% 
  mutate("cluster" = "all_queries")
### add as a separate category "all_queries" to GO_cat_clusters
GO_cat_clusters <- bind_rows(GO_cat_clusters, all_query_genes_cluster) 
tail(GO_cat_clusters)
(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(Gene.gene_name, cluster))

### Whole ubermap as a matrix
#### library genes, e-map score only matrix
ubermap_spread <- ubermap %>% 
  select(Gene_uniq, library.ORF, score) %>% 
  spread(library.ORF, score)
temp_matrix <- as.matrix(ubermap_spread[, -1])
rownames(temp_matrix) <- ubermap_spread$Gene_uniq
query_genes.ubermap.matrix[["all_queries"]] <- temp_matrix


distance_method <- "pearson"
hclust_method <- "complete"
### this function makes subclusters by cutting the tree using cutreeDynamic, 
### it takes the minimal size of the cluster as an argument
## and also deepSplit (from documentation: "The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method")
make_subclusters <- function(input_matrix, GO_cat_name, dendro, dissim, min_size, split_n) { 
  cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, 
                              minClusterSize = min_size, method = "hybrid", deepSplit = split_n, 
                              pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
  clustering.result <- tibble("subcluster" = cut_hybrid, "Gene_uniq" = rownames(input_matrix), "GO" = GO_cat_name)
  unclustered_genes <- clustering.result %>% filter(subcluster == 0)
  clusters <- clustering.result %>% filter(subcluster != 0)
  clusters <- clusters %>% 
    mutate("cluster" = str_c(GO_cat_name, min_size, subcluster, sep = "_")) %>%
    select(cluster, Gene_uniq) %>%
    inner_join(., Gene_uniq_gene_name_index, by = "Gene_uniq") %>%
    arrange(desc(cluster), Gene_uniq)
  return(clusters)
}

### do the clustering for each of the matrices
#pdf(paste0("reverse_corr_of_corr/clusters/", Sys.Date(), "_clusters.pdf"), width = 10)
for (GO_cat in names(query_genes.ubermap.matrix)) {
  temp.mat <- query_genes.ubermap.matrix[[GO_cat]]
  dissim <- Dist(temp.mat, method = distance_method)
  dendro <- hclust(dissim, method = hclust_method)
  #plot(dendro, cex = 0.2, main = paste(GO_cat, distance_method, hclust_method, sep = "_"))
  query_clusters <- tibble("cluster" = character(), "Gene_uniq" = character(), "Gene.gene_name" = character())
  query_clusters <- make_subclusters(input_matrix = temp.mat, GO_cat_name = GO_cat, 
                    dendro = dendro, dissim = dissim, min_size = 30, split_n = 3)
  GO_cat_clusters <- bind_rows(GO_cat_clusters, query_clusters)
}

#dev.off()

(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(Gene_uniq, cluster))

GO_cat_clusters %>% pull(cluster) %>% unique() %>% length()

(cluster_count <- GO_cat_clusters %>% 
    group_by(cluster) %>%
    summarize("count" = n()) %>% 
    arrange(desc(count)))
tail(cluster_count)
cluster_count %>% 
  filter(cluster != "all_queries") %>% 
  ggplot(., aes(x = count)) + geom_histogram() + ggtitle("Distribution of cluster sizes")

##### number of clusters per query gene
(n_clusters_per_protein <- GO_cat_clusters %>% 
  group_by(Gene_uniq, Gene.gene_name) %>%
  summarize("count" = n()) %>%
  ungroup() %>% 
  arrange(desc(count)))
#### how many of the query genes don't belong to any clusters?
n_sorted_query_gene_uniqs <- n_clusters_per_protein %>% pull(Gene_uniq) %>% unique() %>% length()
n_queries <- ubermap %>% pull(Gene_uniq) %>% unique() %>% length()
(unsorted_query_gene_uniqs <- ubermap %>%
  filter(! Gene_uniq %in% n_clusters_per_protein$Gene_uniq) %>% 
  select(Gene_uniq, Gene.gene_name) %>% 
  unique())

## check how many unsorted when I remove the all_queries clusters!
GO_cat_clusters_all_queries_removed <- GO_cat_clusters %>% 
  filter(! grepl("all_queries", cluster))
(n_clusters_per_protein <- GO_cat_clusters_all_queries_removed %>% 
    group_by(Gene_uniq, Gene.gene_name) %>%
    summarize("count" = n()) %>%
    arrange(desc(count)))
(unsorted_query_gene_uniqs <- ubermap %>%
  filter(! Gene_uniq %in% n_clusters_per_protein$Gene_uniq) %>% 
  select(Gene_uniq, Gene.gene_name) %>% 
  unique())
### if I don't consider the all_queries clusters there are 564 unsorted genes
### many of them seem to be of either unknown or poorly defined function

outfilename <- str_c("reverse_corr_of_corr/clusters/GO_slims", Sys.Date(), "pearson_complete_query_genes_clusters.txt", sep = "_")
write_tsv(GO_cat_clusters, path = outfilename)


