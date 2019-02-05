library(tidyverse)
library(factoextra)
library(fastcluster)
data <- read_tsv("basic_E-MAP_data/20180921_preprocessed_ubermap_500_overlap_all.txt")

orfs_to_keep <- c("YLR293C", "YGL097W", "YMR235C", "YDR002W", "YJR074W", "YGR218W", "YLR347C", "YDR335W")
high_score_lib_genes <- read_tsv("library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt")

ubermap <- data %>% 
  filter(ORF %in% orfs_to_keep) %>% 
  filter(! is.na(score)) %>% 
  mutate("Gene.gene_name" = ifelse(Gene.gene_name == "GSP1", Gene_uniq, Gene.gene_name))
ubermap %>% 
  filter(ORF != "YLR293C") %>% 
  group_by(library.ORF) %>% 
  summarise("count" = n()) %>% 
  ggplot(aes(x = count)) + geom_histogram()
# lib.genes.outside.box <- ubermap %>%
#   group_by(library.ORF) %>% 
#   summarize("abs.max" = max(abs(score), na.rm = T)) %>% 
#   arrange(abs.max) %>% 
#   filter(abs.max > 2) %>% 
#   pull(library.ORF)
library.ORFs_with_overlap <- ubermap %>% 
  filter(ORF != "YLR293C") %>% 
  #filter(library.ORF %in% lib.genes.outside.box) %>% 
  group_by(library.ORF) %>% 
  summarise("count" = n()) %>% 
  filter(count >= 3) %>% 
  pull(library.ORF) %>% unique()

ubermap <- ubermap %>% 
  filter(library.ORF %in% library.ORFs_with_overlap)

ubermap_spread <- ubermap %>% 
  mutate('residue' = as.numeric(substring(Gene.gene_name, 2, (nchar(Gene.gene_name)-1)))) %>% 
  select(residue, Gene.gene_name, library.gene_name, score) %>% 
  spread(library.gene_name, score)
ubermap_spread <- Filter(function (x) !all (is.na(x)), ubermap_spread)  # this removes columns that are all NA
write_tsv(ubermap_spread, path = "clustering_muts_and_core_regulators/emap_with_core_regs.txt")
ubermap_spread <- data.frame(ubermap_spread)
rownames(ubermap_spread) <- ubermap_spread$Gene.gene_name
ubermap_spread_matrix <- ubermap_spread[, -c(1,2)]
ubermap_mutants_hclust <- hcut(dist(ubermap_spread_matrix), k = 6)
fviz_dend(ubermap_mutants_hclust, main = "mutants and core regulators clustered by E-MAP score k = 6", cex = 0.6)  ## from ‘factoextra’ package 
### order of mutants/partners based on hierarchical clustering
order_mutants <- ubermap_mutants_hclust$labels[ubermap_mutants_hclust$order]
order_residues <- ubermap_spread$residue[ubermap_mutants_hclust$order]
## interfaces
interfaces <- read_tsv("SASA_interfaces.txt")
interface_residues_completed <- interfaces %>%
  filter(interface == "core") %>% 
  select(partner, yeastresnum, deltarASA) %>%
  complete(partner, nesting(yeastresnum), fill = list(deltarASA = 0)) %>% 
  arrange(yeastresnum) %>% 
  filter(yeastresnum %in% order_residues & partner %in% order_mutants)
partners_add_on <- interface_residues_completed 



orfs_to_keep <- c("YLR293C")
ubermap <- data %>% 
  filter(ORF %in% orfs_to_keep) %>% 
  filter(! is.na(score)) %>% 
  mutate("Gene.gene_name" = ifelse(Gene.gene_name == "GSP1", Gene_uniq, Gene.gene_name))
ubermap_spread <- ubermap %>% 
  select(Gene.gene_name, library.gene_name, score) %>% 
  spread(library.gene_name, score)
ubermap_spread <- Filter(function (x) !all (is.na(x)), ubermap_spread)  # this removes columns that are all NA
write_tsv(ubermap_spread, path = "clustering_muts_and_core_regulators/emap_with_core_regs.txt")
ubermap_spread <- data.frame(ubermap_spread)
rownames(ubermap_spread) <- ubermap_spread$Gene.gene_name
ubermap_spread <- ubermap_spread[, -1]
ubermap_mutants_hclust <- hcut(dist(ubermap_spread), k = 6)
fviz_dend(ubermap_mutants_hclust, main = "mutants and core regulators clustered by E-MAP score k = 6", cex = 0.6)  ## from ‘factoextra’ package 
### order of mutants based on hierarchical clustering
ubermap_mutants_hclust$labels[ubermap_mutants_hclust$order]









### correlation matrix for core subset
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}





