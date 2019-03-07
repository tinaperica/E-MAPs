library(tidyverse)
set.seed(2019)
### this function is defined from the original GSEA paper R code (Subramanian, A. et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102, 15545–15550 (2005).)
GSEA.EnrichmentScore <- function(gene.list, gene.set, weight = 1, correl.vector = NULL) {
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  alpha <- weight
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))  
}
  
Fast.GSEA.EnrichmentScore <- function(gene.list, gene.set, weight = 1, correl.vector = NULL) {  
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  loc.vector <- vector(length = N, mode = "numeric")
  peak.res.vector <- vector(length = Nh, mode = "numeric")
  valley.res.vector <- vector(length = Nh, mode = "numeric")
  tag.correl.vector <- vector(length = Nh, mode = "numeric")
  tag.diff.vector <- vector(length = Nh, mode = "numeric")
  tag.loc.vector <- vector(length = Nh, mode = "numeric")
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  if (weight == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weight == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weight == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weight
    tag.correl.vector <- abs(tag.correl.vector)
  }
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits = 5)
  return(list(ES = ES))
}

kendall_corr <- read_tsv("integrating_biophysical_parameters/kendall_correlation.txt") %>% 
  select(gene_name, "GAP_kcat_Km", "GEF_kcat_Km")
library_gene_names <- kendall_corr %>% pull(gene_name) %>% unique()

gene_groups <- read_tsv("GSEA_like_analysis/gene_groups.txt") %>% 
  filter(Gene %in% library_gene_names)
larger_than_3_terms <- gene_groups %>%    #### still larger than 3 after overlapping with emap data
  group_by(term) %>% summarise("count" = n()) %>% 
  filter(count > 3) %>% pull(term) %>% unique()
gene_groups <- gene_groups %>% 
  filter(term %in% larger_than_3_terms) %>% 
  select("gene_name" = Gene, term)


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
get_reordered_cormat <- function(data) {
  cormat <- round(cor(data[, -1], use = "pairwise.complete.obs", method = "spearman"), 2)
  cormat <- reorder_cormat(cormat)
  return(cormat)
}
get_order_by_corr <- function(data) {
  cormat <- get_reordered_cormat(data)
  ordered <- rownames(cormat)
}

ordered_genes <- kendall_corr %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  gather(param, value, -gene_name) %>% 
  spread(gene_name, value) %>% 
  get_order_by_corr()

kendall_corr %>% 
  select(gene_name, GAP_kcat_Km, GEF_kcat_Km) %>% 
  mutate("gene_name" = factor(gene_name, ordered_genes)) %>% 
  rename("GAP" = GAP_kcat_Km, "GEF" = GEF_kcat_Km) %>% 
  gather(param, value, -gene_name) %>% 
  ggplot(aes(param, gene_name, fill = value)) + geom_tile() +
  scale_fill_viridis() + 
  xlab("Kendall correlation between genetic interactions and GTPase cycle") +
  ylab("S. cerevisae gene") + theme(axis.text.y = element_blank())

