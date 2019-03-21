### this script makes hierarchical clusters of query genes based on GO slims terms as well
# as on correlations of genetic interactions (ubermap e-map scores)

library(fastcluster)  # overwrites the hclust function and makes it faster
library(dynamicTreeCut)   ### decides where to cut the hclust tree to make clusters
library(amap)
library(tidyverse)

### ORF to gene name annotation from SGD
orf_gene_name_index <- read_delim("orf_gene_GO_sgd_annotation.txt", col_names = F, delim = "\t")
orf_index <- unique(tibble("ORF" = orf_gene_name_index$X1, "gene_name" = orf_gene_name_index$X2))
rm(orf_gene_name_index)
#############################
##### read the matcher of unique query identifiers with gene names/ ORFs
matcher <- read_tsv("basic_E-MAP_data/20190101_SGAubermap_gene_uniq_matcher.txt")

###### load SGAubermap correlations data - loads an SGA_correlations tibble
load("corr_of_corr_with_SGA/20190102_all_correlations.RData")
### add the ORF columns
SGA_correlations <- matcher %>%
  select(query_uniq, query_ORF) %>% 
  inner_join(., SGA_correlations, by = c("query_uniq" = "query_uniq1")) %>% 
  rename("query_uniq1" = query_uniq, "query_ORF1" = query_ORF) %>% 
  filter(query_ORF1 != "YLR293C")
query_uniq_ORF_index <- SGA_correlations %>% 
  select(query_uniq1, query_ORF1) %>% 
  unique()
### first make clusters based on GO slims
GO_slims <- read_delim( "reverse_corr_of_corr/20180921_go_slim_mapping.tab.txt", col_names = F, delim = "\t") %>% 
  select("ORF" = X1, "Gene" = X2, "GO_Slim_term" = X5) %>% 
  mutate("Gene" = ifelse(is.na(Gene), ORF, Gene))
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
### by pearson correlation of CORRELATIONS
### also have clusters where there was no previous subclustering of GO slims categories 


query_genes.SGAcorr.matrix <- list()  # all the matrices are stored in a single list
GO_cat_clusters <- tibble("cluster" = character(), "query_uniq1" = character(), "query_ORF1" = character())
for (cat in names(GO_slims_for_clusters)) {
  GO_slims_cat <- GO_slims_for_clusters[[cat]]
  GO_slims_ORFs <- GO_slims %>% 
    filter(GO_Slim_term %in% GO_slims_cat) %>%
    pull("ORF") %>% unique()
  clust.SGAcorr_queries <- SGA_correlations %>% 
    filter(query_ORF1 %in% GO_slims_ORFs)
  clust_unique_queries <- clust.SGAcorr_queries %>% 
    pull(query_uniq1) %>% unique() %>% sort()
  if (length(clust_unique_queries) >= 50) {
    cluster <- clust.SGAcorr_queries %>% 
      select(query_uniq1, query_ORF1) %>% 
      mutate("cluster" = cat) %>% unique() %>% 
      select(cluster, query_uniq1, query_ORF1) %>% 
      arrange(cluster, query_uniq1)
    GO_cat_clusters <- bind_rows(GO_cat_clusters, cluster)
    clust_SGAcorr_spread <- clust.SGAcorr_queries %>% 
      select(query_uniq1, query_uniq2, w_correlation) %>% 
      spread(query_uniq2, w_correlation)
    temp_matrix <- as.matrix(clust_SGAcorr_spread[, -1])
    rownames(temp_matrix) <- clust_SGAcorr_spread$query_uniq1
    ### add each of these clusters onto the query_genes.ubermap.matrix so I can subscluster them with hclust later
    query_genes.SGAcorr.matrix[[str_c(cat, "GO", sep = "_")]] <- temp_matrix
  }
}

# check that it looks as expected
head(query_genes.SGAcorr.matrix[["ribosome_GO"]])[,1:20]

outfilename <- str_c("corr_of_corr_with_SGA/clusters/broad_GO_cat_clusters_", Sys.Date(), ".txt", sep = "")
write_tsv(GO_cat_clusters, path = outfilename)

### add the total of all the query genes as a "cluster"
all_query_genes_cluster <- SGA_correlations %>% 
  select(query_uniq1, query_ORF1) %>% 
  unique() %>% 
  mutate("cluster" = "all_queries")
### add as a separate category "all_queries" to GO_cat_clusters
GO_cat_clusters <- bind_rows(GO_cat_clusters, all_query_genes_cluster) 
tail(GO_cat_clusters)
(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(query_uniq1, cluster))

### Whole SGA_correlations as a matrix
SGA_corr_spread <- SGA_correlations %>% 
  select(query_uniq1, query_uniq2, w_correlation) %>% 
  spread(query_uniq2, w_correlation)
temp_matrix <- as.matrix(SGA_corr_spread[, -1])
rownames(temp_matrix) <- SGA_corr_spread$query_uniq1
query_genes.SGAcorr.matrix[["all_queries"]] <- temp_matrix


distance_method <- "pearson"
hclust_method <- "complete"
### this function makes subclusters by cutting the tree using cutreeDynamic, 
### it takes the minimal size of the cluster as an argument
## and also deepSplit (from documentation: "The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method")
make_subclusters <- function(input_matrix, GO_cat_name, dendro, dissim, min_size, split_n) { 
  cut_hybrid <- cutreeDynamic(dendro = dendro, cutHeight = NULL, 
                              minClusterSize = min_size, method = "hybrid", deepSplit = split_n, 
                              pamStage = T, distM = as.matrix(dissim), maxPamDist = 0, verbose = 0)
  clustering.result <- tibble("subcluster" = cut_hybrid, "query_uniq1" = rownames(input_matrix), "GO" = GO_cat_name)
  unclustered_genes <- clustering.result %>% filter(subcluster == 0)
  clusters <- clustering.result %>% filter(subcluster != 0)
  clusters <- clusters %>% 
    mutate("cluster" = str_c(GO_cat_name, min_size, subcluster, sep = "_")) %>%
    select(cluster, query_uniq1) %>%
    inner_join(., query_uniq_ORF_index, by = "query_uniq1") %>%
    arrange(desc(cluster), query_uniq1)
  return(clusters)
}

### do the clustering for each of the matrices
#pdf(paste0("reverse_corr_of_corr/clusters/", Sys.Date(), "_clusters.pdf"), width = 10)
for (GO_cat in names(query_genes.SGAcorr.matrix)) {
  temp.mat <- query_genes.SGAcorr.matrix[[GO_cat]]
  dissim <- Dist(temp.mat, method = distance_method)
  dendro <- hclust(dissim, method = hclust_method)
  query_clusters <- tibble("cluster" = character(), "query_uniq1" = character(), "query_ORF1" = character())
  query_clusters <- make_subclusters(input_matrix = temp.mat, GO_cat_name = GO_cat, 
                    dendro = dendro, dissim = dissim, min_size = 30, split_n = 3)
  GO_cat_clusters <- bind_rows(GO_cat_clusters, query_clusters)
}

#dev.off()

(GO_cat_clusters <- GO_cat_clusters %>% 
    arrange(query_uniq1, cluster))

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
  group_by(query_uniq1) %>%
  summarize("count" = n()) %>%
  ungroup() %>% 
  arrange(desc(count)))
#### how many of the query genes don't belong to any clusters?
n_sorted_query_gene_uniqs <- n_clusters_per_protein %>% pull(query_uniq1) %>% unique() %>% length()
n_queries <- SGA_correlations %>% pull(query_uniq1) %>% unique() %>% length()
(unsorted_query_gene_uniqs <- SGA_correlations %>%
  filter(! query_uniq1 %in% n_clusters_per_protein$query_uniq1) %>% 
  select(query_uniq1) %>% 
  unique())

## check how many unsorted when I remove the all_queries clusters!
GO_cat_clusters_all_queries_removed <- GO_cat_clusters %>% 
  filter(! grepl("all_queries", cluster))
(n_clusters_per_protein <- GO_cat_clusters_all_queries_removed %>% 
    group_by(query_uniq1) %>%
    summarize("count" = n()) %>%
    arrange(desc(count)))
(unsorted_query_gene_uniqs <- SGA_correlations %>%
  filter(! query_uniq1 %in% n_clusters_per_protein$query_uniq1) %>% 
  select(query_uniq1) %>% 
  unique())
### if I don't consider the all_queries clusters there are 685 unsorted genes
### many of them seem to be of either unknown or poorly defined function

outfilename <- str_c("corr_of_corr_with_SGA/clusters/GO_slims", Sys.Date(), "pearson_complete_query_genes_clusters.txt", sep = "_")
write_tsv(GO_cat_clusters, path = outfilename)


