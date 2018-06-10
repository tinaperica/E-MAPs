# get a distribution of number of genes in a cluster
# sample random numbers from that distribution and make random clusters of that size
# make 45 such random clusters, 3 independent times
### check the distribution of number of clusters each gene is
### if that distribution is not similar to the distribution for the real clusters, look for cluster hub proteins
## i.e proteins that are in many clusters and give them some kind of higher weight, so that they appear in multiple clusters in random data too
### then run the full corr of corr analysis with the 3 random cluster sets and compare the heatmaps with real data
library(tidyverse)
library(RColorBrewer)
number_of_simulations <- 3
#### these are the current clusters I am working with
clusters.data <- read_delim(file = "2018-01-17_GO_slims_pearson_complete_clusters.txt", delim = "\t")

(genes.per.cluster <- clusters.data %>%
  group_by(cluster) %>%
  summarise("count" = n()) )

(clusters.per.gene <- clusters.data %>%
  group_by(gene_name) %>%
  summarise("count" = n()) )

### define some basic numbers to use throughout
###### simulated data will keep the same number of clusters (and keep the names)
clusters <- unique(clusters.data$cluster)
number_of_clusters <- length(clusters)
genes <- unique(clusters.data$gene_name)
orfs <- unique(clusters.data$ORF)
number_of_genes <- length(genes)
gene_to_orf_index <- clusters.data %>% 
        select(c("gene_name", "ORF")) %>%
        group_by(gene_name, ORF) %>%
        filter(row_number() == 1)

### this part simulates the sizes of clusters and the number of clusters each gene is in
simulated_cluster_sizes <- list()
cluster.density.model <- density(genes.per.cluster$count)
gene.density.model <- density(clusters.per.gene$count, n = 2^10)
colors <- brewer.pal(n = number_of_simulations, "Set1")
plot(cluster.density.model, main = "Distribution of cluster sizes", ylim = c(0, 0.025))
for (i in 1:number_of_simulations) {
  set.seed(1984 + i)
  simulated_cluster_sizes[["cluster_sizes"]][[i]] <- tibble("cluster" = clusters, 
          "size" = round(sample(cluster.density.model$x, 
          size = number_of_clusters, replace = T, prob = cluster.density.model$y)))
  lines(density(simulated_cluster_sizes[["cluster_sizes"]][[i]][["size"]]), col = colors[i])
  legend("topright", fill = colors, legend = c("simulation 1", "simulation 2", "simulation 3"))
}
plot(gene.density.model, main = "Distribution of number of clusters per gene", ylim = c(0, 0.65))
for (i in 1:number_of_simulations) {
  set.seed(1984 + i)
  simulated_cluster_sizes[["clusters_per_gene"]][[i]] <- tibble("orf" = orfs, 
          "n_clusters" = round(sample(gene.density.model$x, 
          size = number_of_genes, replace = T, prob = gene.density.model$y)))
  lines(density(simulated_cluster_sizes[["clusters_per_gene"]][[i]][["n_clusters"]]), col = colors[i])
  legend("topright", fill = colors, legend = c("simulation 1", "simulation 2", "simulation 3"))
}
#######################
#pdf(file = "clustered_correlations/clusters/simulated_distributions.pdf")
for (sim_n in 1:number_of_simulations) {
  simulated_clusters <- tibble(ORF = character(), cluster = character(), gene_name = character())
  filename <- str_c("simulated_cluster_", sim_n, ".txt")
  sim_n_clusters <- simulated_cluster_sizes[["cluster_sizes"]][[sim_n]]
  sim_n_orfs <- simulated_cluster_sizes[["clusters_per_gene"]][[sim_n]]
  ##### populate one by one cluster
  sim_n_genes_to_sample_from <- sim_n_orfs %>% filter(n_clusters > 0)
  for (cn in 1:number_of_clusters) {
    clust <- sim_n_clusters$cluster[cn]
    clust_size <- sim_n_clusters[["size"]][sim_n_clusters$cluster == clust]
    ##### sim_n_genes_to_sample runs out of geens before all the clusters are populated
    #### to hack around that problem, once the pool of genes in sim_n_genes_to_sample becomes too small, add some more random genes
    if (nrow(sim_n_genes_to_sample_from) < clust_size) {
      diff <- clust_size - nrow(sim_n_genes_to_sample_from)
      sim_n_genes_to_sample_from <- add_row(sim_n_genes_to_sample_from, orf = sample(sim_n_orfs$orf, size = diff), n_clusters = 1)
    }
    temp_orfs <- sample(sim_n_genes_to_sample_from$orf, size = sim_n_clusters$size[sim_n_clusters$cluster == clust], replace = F)
    temp_tibble <- tibble(cluster = clust, ORF = temp_orfs)
    temp_tibble <- merge(temp_tibble, gene_to_orf_index)
    simulated_clusters <- add_row(simulated_clusters, cluster = temp_tibble$cluster, ORF = temp_tibble$ORF, 
                          gene_name = temp_tibble$gene_name)
    sim_n_genes_to_sample_from <- sim_n_genes_to_sample_from %>% 
                  mutate(n_clusters = if_else(orf %in% temp_orfs, n_clusters - 1, n_clusters)) %>%
                  filter(n_clusters > 0)
  }
  
  sim.genes.per.cluster <- simulated_clusters %>% group_by(cluster) %>% summarise("count" = n())
  sim.genes.per.cluster.density.model <- density(sim.genes.per.cluster$count)
  plot(cluster.density.model, main = str_c("Distribution of cluster sizes ", sim_n), ylim = c(0, 0.025), col = colors[1])
  legend("topright", fill = colors, legend = c("real clusters", "simulated distribution", "simulated data"))
  lines(density(simulated_cluster_sizes[["cluster_sizes"]][[sim_n]][["size"]]), col = colors[2], lwd = 2)
  lines(sim.genes.per.cluster.density.model, col = colors[3], lty = 2)
  sim.clusters.per.gene <- simulated_clusters %>% group_by(gene_name) %>% summarise("count" = n())
  sim.clusters.per.gene.density.model <- density(sim.clusters.per.gene$count)
  plot(gene.density.model, main = str_c("Distribution of number of clusters per gene ", sim_n), col = colors[1])
  legend("topright", fill = colors, legend = c("real clusters", "simulated distribution", "simulated data"))
  lines(density(simulated_cluster_sizes[["clusters_per_gene"]][[sim_n]][["n_clusters"]]), col = colors[2], lwd = 2)
  lines(sim.clusters.per.gene.density.model, col = colors[3], lty = 2)
  write.table(simulated_clusters, file = filename, quote = F, row.names = F, sep = "\t")
}
#dev.off()




