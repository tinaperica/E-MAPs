library(tidyverse)
library(factoextra)
library(fastcluster)
library(RColorBrewer)
library(ggrepel)

date <- as.character(Sys.Date())
### load in e.map data
load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
e.map <- as_tibble(e.map) %>% 
  mutate("mutant" = as.character(mutant), library_ORF = as.character(library_ORF), 
         "library_gene_name" = as.character(library_gene_name)) %>% 
  mutate("mutant" = ifelse(mutant == "GSP1-NAT", "WT", mutant)) 
#### simple hierarchical clustering of mutants based on E-MAP scores

### input all the parameters for mutants
GEF_kin <- read_tsv("hierarchical_clustering_of_parameters/GEF_kinetics_param.txt", col_names = T) %>% 
  gather("measure", "value", -mutant)
GAP_kin <- read_tsv("hierarchical_clustering_of_parameters/GAP_MM_parameters_clean.txt", col_names = T) %>% 
  gather("measure", "value", -mutant)
NMR_data <- read_tsv("hierarchical_clustering_of_parameters/20190224_nmr_data.txt", col_names = T) %>% 
  select("mutant" = "gsp1_mutation", "gamma2" = gamma2_percent) %>% 
  gather("measure", "value", -mutant)
apms_foldchange <- read_tsv("hierarchical_clustering_of_parameters/apms_log2_fold_change.txt", col_names = T) %>% 
  select(mutant, log2FC, tag, gene_name) %>% 
  filter(tag == "N") %>% 
  group_by(mutant, gene_name) %>% 
  mutate("avg_log2FC" = mean(log2FC)) %>% 
  ungroup() %>% 
  select(mutant, "measure" = gene_name, "value" = avg_log2FC) %>% 
  filter(measure %in% c("SRM1", "RNA1", "YRB1", "MOG1"))
all_parameters <- bind_rows(GEF_kin, GAP_kin, NMR_data, apms_foldchange)
all_parameters <- all_parameters %>% 
  ##### add log2fold change of zero for WT
  bind_rows(., tibble("measure" = c("SRM1", "RNA1", "YRB1", "MOG1"), mutant = "WT", value = rep(0, 4))) %>% 
  unique() %>% 
  group_by(measure) %>% 
  mutate("raw_value" = value) %>% 
  mutate("value" = scale(raw_value)) %>% 
  ungroup()
###########


### hierarchical clustering of the emap data
spread_emap <- e.map %>% 
  filter(mutant %in% unique(all_parameters$mutant)) %>% 
  select(library_gene_name, mutant, score) %>% 
  spread(library_gene_name, score)
spread_emap <- data.frame(spread_emap)
rownames(spread_emap) <- spread_emap$mutant
spread_emap <- spread_emap[, -1]
spread_emap <- Filter(function (x) !all (is.na(x)), spread_emap)  # this removes columns that are all NA
head(spread_emap[1:10, 1:10])
hclust_mut_data <- hcut(spread_emap, k = 4)
order_mutants <- hclust_mut_data$labels[hclust_mut_data$order]
fviz_dend(hclust_mut_data, k_colors = brewer.pal(4, "Set1"),
          main = "mutants clustered by E-MAP score, k = 4")  # from ‘factoextra’ package 




cluster_by_parameters <- function(emap, biophy_params, k) {
  plots <- list()
  param_subset <- all_parameters %>% 
    filter(measure %in% biophy_params)
  mutants_with_all_params <- param_subset %>% 
    select(mutant, measure) %>% 
    group_by(mutant) %>% 
    summarise("n_param" = n()) %>% 
    ungroup() %>% 
    mutate("max_n_param" = max(n_param)) %>% 
    filter(n_param == max_n_param) %>% 
    pull(mutant)
  emap.subset <- emap %>% 
    filter(mutant %in% mutants_with_all_params)
  
  ### hierarchical clustering of the emap data
  spread_emap <- emap.subset %>% 
    #filter(mutant %in% unique(all_parameters$mutant)) %>% 
    select(library_gene_name, mutant, score) %>% 
    spread(library_gene_name, score)
  spread_emap <- data.frame(spread_emap)
  rownames(spread_emap) <- spread_emap$mutant
  spread_emap <- spread_emap[, -1]
  spread_emap <- Filter(function (x) !all (is.na(x)), spread_emap)  # this removes columns that are all NA
  head(spread_emap[1:10, 1:10])
  hclust_emap_data <- hcut(spread_emap, k = k)
  
  mut_groups <- data.frame(cutree(hclust_emap_data, k = k))
  names(mut_groups) <- "group"
  mut_groups <- tibble("mutant" = rownames(mut_groups), "group" = mut_groups$group) %>% 
    inner_join(., tibble("group" = seq(1, k, 1), "color" = brewer.pal(k, "Set1")), by = "group") %>% 
    arrange(group)
  order_mutants <- hclust_emap_data$labels[hclust_emap_data$order]
  mut_groups <- mut_groups %>% 
    mutate("mut_fact" = factor(mutant, order_mutants)) %>% 
    arrange(mut_fact)
  plots[["by_emap"]] <- fviz_dend(hclust_emap_data, label_cols =  mut_groups$color, k_colors = "black",
            main = str_c("mutants clustered by E-MAP score, k = ", k))  # from ‘factoextra’ package 
  
  #### cluster by parameters
  param_subset_spread <- param_subset %>% 
    filter(mutant %in% mutants_with_all_params) %>% 
    select(-raw_value) %>% 
    spread(measure, value)
  param_subset_spread <- data.frame(param_subset_spread)
  rownames(param_subset_spread) <- param_subset_spread$mutant
  param_subset_spread <- param_subset_spread[, -1]
  hclust_param_data <- hcut(param_subset_spread, k = k)
  order_mutants_by_param <- hclust_param_data$labels[hclust_param_data$order]
  mut_groups <- mut_groups %>% 
    mutate("mut_fact" = factor(mutant, order_mutants_by_param)) %>% 
    arrange(mut_fact)
  plots[["by_param"]] <- fviz_dend(hclust_param_data, label_cols =  mut_groups$color, k_colors = "black", 
            main = str_c( c("mutants clustered by: (", biophy_params, ") -> k =", k), collapse = " "))
  pdf(file = str_c(c("hierarchical_clustering_of_parameters/", date, "_", biophy_params, ".pdf"), collapse = ""), width = 10)
  print(plots)
  dev.off()
}

available_parameters <- all_parameters %>% pull(measure) %>% unique()
cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km"), k = 6)
cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat_Km", "GEF_kcat_Km"), k = 4)

cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2"), k = 4)

cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1"), k = 4)
cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "MOG1"), k = 4)
cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "YRB1"), k = 4)
cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "MOG1", "YRB1"), k = 4)
cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2", "SRM1", "RNA1"), k = 4)
cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2", "SRM1", "RNA1", "MOG1", "YRB1"), k = 4)



WT_parameters <- all_parameters %>% 
  filter(mutant == "WT") %>% 
  select(-value, -mutant) %>% 
  rename("wt_value" = raw_value)
all_parameters %>% filter(measure %in% c("gamma2", "GAP_kcat_Km")) %>% 
  select(-value) %>% 
  inner_join(., WT_parameters, by = c("measure")) %>% 
  mutate("rel_value" = raw_value/wt_value) %>% 
  select(-raw_value, -wt_value) %>% 
  spread(measure, rel_value) %>% 
  ggplot(aes(x = gamma2, y = log(GAP_kcat_Km), label = mutant)) +
  geom_point() +
  geom_text_repel() +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  xlab("% state 2") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave(filename = "hierarchical_clustering_of_parameters/GAP_enzyme_eff_vs_gamma2.pdf", width = 12)


rel_GAP_GEF_efficiency <- all_parameters %>% 
  filter(measure %in% c("GEF_kcat_Km", "GAP_kcat_Km")) %>% 
  select(-value) %>% 
  inner_join(., WT_parameters, by = c("measure")) %>% 
  mutate("rel_value" = raw_value/wt_value) %>% 
  select(-raw_value, -wt_value) 
  
  
rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  ggplot(aes(x = log(GEF_kcat_Km), y = log(GAP_kcat_Km), label = mutant)) +
  geom_point() +
  geom_text_repel() +
  xlab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GEF mediated nucleotide exchange") +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

ggsave(filename = "hierarchical_clustering_of_parameters/GAP_enzyme_eff_vs_GEF_enzyme_eff.pdf", width = 12)

GAP_GEF_ratio <- rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  mutate("GAP_GEF_ratio" = GAP_kcat_Km / GEF_kcat_Km) %>% 
  filter(! is.na(GAP_GEF_ratio)) %>% 
  arrange(GAP_GEF_ratio) 
mutants_by_GAP_GEF_ratio <- GAP_GEF_ratio %>% 
  pull(mutant) %>% unique()

GAP_GEF_ratio %>% 
  mutate("mutant" = factor(mutant, mutants_by_GAP_GEF_ratio)) %>% 
  ggplot(aes(x = mutant, y = log(GAP_GEF_ratio))) + 
  geom_point(size = 6, alpha = 0.5) +
  ylab("ln(relative GAP over relative GEF efficiency)") +
  xlab("Gsp1 point mutant") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
        axis.title = element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", alpha = 0.5)



# 
# GEF_kin <- data.frame(GEF_kin)
# rownames(GEF_kin) <- GEF_kin$mutant
# GEF_kin <- GEF_kin[, -1]
# k = 8
# hclust_kin_data <- hcut(GEF_kin, k = k)
# # fviz_dend(hclust_kin_data,  k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), 
# #           main = "mutants clustered by nucleotide exchange kinetics, k = 4")
# fviz_dend(hclust_kin_data,  k_colors = brewer.pal(k, "Set1"), 
#       main = str_c("mutants clustered by nucleotide exchange kinetics, k = ", k))
# GEF_groups <- data.frame(cutree(hclust_kin_data, k = k))
# names(GEF_groups) <- "group"
# GEF_groups <- tibble("mutant" = rownames(GEF_groups), "group" = GEF_groups$group) %>% 
#   inner_join(., tibble("group" = seq(1, k, 1), "color" = brewer.pal(k, "Set2")), by = "group") %>% 
#   arrange(group)
# order_mutants <- hclust_kin_data$labels[hclust_kin_data$order]
# GEF_groups <- GEF_groups %>% 
#   mutate("mut_fact" = factor(mutant, order_mutants)) %>% 
#   arrange(mut_fact)
# fviz_dend(hclust_kin_data,  label_cols =  GEF_groups$color, k_colors = "black", 
#           main = str_c("mutants clustered by nucleotide exchange kinetics, k = ", k))
# spread_emap <- e.map %>% 
#   filter(mutant %in% rownames(GEF_kin)) %>% 
#   select(library_gene_name, mutant, score) %>% 
#   spread(library_gene_name, score)
# spread_emap <- data.frame(spread_emap)
# rownames(spread_emap) <- spread_emap$mutant
# spread_emap <- spread_emap[, -1]
# spread_emap <- Filter(function (x) !all (is.na(x)), spread_emap)  # this removes columns that are all NA
# head(spread_emap[1:10, 1:10])
# hclust_mut_data <- hcut(spread_emap, k = 4)
# order_mutants <- hclust_mut_data$labels[hclust_mut_data$order]
# GEF_groups <- GEF_groups %>% 
#   mutate("mut_fact" = factor(mutant, order_mutants)) %>% 
#   arrange(mut_fact)
# fviz_dend(hclust_mut_data, label_cols =  GEF_groups$color, k_colors = "black",
#           main = "mutants clustered by E-MAP score, clusters by GEF kinetic parameters")  # from ‘factoextra’ package 
# 
# 
