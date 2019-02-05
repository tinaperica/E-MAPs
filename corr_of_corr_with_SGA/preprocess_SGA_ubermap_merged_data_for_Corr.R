library(tidyverse)
library(ggrepel)
library(ggcorrplot)
library(gmp)
library(GGally)
library(reshape2)
abs.pmax <- function(x, y) {
  df <- data.frame(x, y)
  pmax <- vector()
  for (i in 1:length(x)) {
    pair <- df[i, c("x", "y")]
    max <- pair[ which.max(abs(pair)) ][1,1]
    pmax <- append(pmax, max)
  }
  return(pmax)
}
asym_pnorm_cropped <- function(data) {  ### return 1 - 2*probability that the value is same or less (for neg scores) or same or greater for pos scores
  prob_weights <- vector()
  for (x in data) {
    if (x < -3) {
       p <- pnorm(x, mean = 0, sd = 5)
    } else {
      if (x > 2) {
        p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      } else {
        p <- 0.5
      }
    }
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}
asym_pnorm <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
    } 
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}
asym_pnorm_signed <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
      prob_weights <- append(prob_weights, -1 * (1 - 2*p))
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      prob_weights <- append(prob_weights, 1 - 2*p)
    } 
  }
  return(prob_weights)
}
asym_pnorm_signed_cropped <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < -3) {
      p <- pnorm(x, mean = 0, sd = 5)
      prob_weights <- append(prob_weights, -1 * (1 - 2*p))
    } else {
      if (x > 2) {
        p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
        prob_weights <- append(prob_weights, 1 - 2*p)
      } else {
        prob_weights <- append(prob_weights, 0)
      }
    } 
  }
  return(prob_weights)
}
w_pearson <- function( x, y) {
  ## pick pairwise observable
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y) # max absolute value E-MAP score of the two
  w_cropped <- asym_pnorm_cropped(max)  ### 0 for E-MAP scores around 0, used for filtering
  w <- asym_pnorm(max)   ## get the weight for the pair
  df <- cbind(df, data.frame(w, w_cropped))
  df_clean <- df[df$w_cropped != 0, ]  ## df_clean is only for filtering
  #if (length(df_clean$x) > 1 & length(df$x) > 3) {
  if (length(df_clean$x) > 3 & length(df$x) > 5) {
    # x and y weighted means
    mean_x <- sum(x * w) / sum(w)
    mean_y <- sum(y * w) / sum(w)
    # Compute the weighted variance
    vx <- sum( w * (x - mean_x)^2 ) / sum(w)
    vy <- sum( w * (y - mean_y)^2 ) / sum(w)
    # Compute the covariance
    vxy <- sum( w * (x - mean_x) * (y - mean_y) ) / sum(w)
    # Compute the correlation
    w_correlation <- (vxy / sqrt(vx * vy))
  } else {
    w_correlation <- NaN
  }
  return(w_correlation)
}
w_pearson_cropped <- function( x, y) {
  ## pick pairwise observable
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y) # max absolute value E-MAP score of the two
  w <- asym_pnorm_cropped(max)  ### 0 for E-MAP scores around 0, used for filtering
  if (length(w[w != 0]) > 5) {
    # x and y weighted means
    mean_x <- sum(x * w) / sum(w)
    mean_y <- sum(y * w) / sum(w)
    # Compute the weighted variance
    vx <- sum( w * (x - mean_x)^2 ) / sum(w)
    vy <- sum( w * (y - mean_y)^2 ) / sum(w)
    # Compute the covariance
    vxy <- sum( w * (x - mean_x) * (y - mean_y) ) / sum(w)
    # Compute the correlation
    w_correlation <- (vxy / sqrt(vx * vy))
  } else {
    w_correlation <- NaN
  }
  return(w_correlation)
}
SoftCosSim <- function(x, y) {
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  df$wx <- asym_pnorm_signed(df$x)
  df$wy <- asym_pnorm_signed(df$y)
  max <- abs.pmax(df$x, df$y)
  w_cropped <- asym_pnorm_cropped(max)
  df <- cbind(df, data.frame(w_cropped))
  df_clean <- df[df$w_cropped != 0, ]
  df$w_sim <- (2 - abs(df$wx - df$wy)) / 2
  x <- df$x
  y <- df$y
  w_sim <- df$w_sim
  if (length(df_clean$x) > 3 & length(df$x) > 5) {
    soft_cos_sim <-   sum(w_sim * x * y) / ( sqrt( sum(w_sim * x^2) ) * sqrt( sum(w_sim * y^2) ) )
  } else {
    soft_cos_sim <- NaN
  }
  return(soft_cos_sim)
}

SoftCosSim_cropped <- function(x, y) {
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  df$wx <- asym_pnorm_signed_cropped(df$x)
  df$wy <- asym_pnorm_signed_cropped(df$y)
  max <- abs.pmax(df$x, df$y)
  #w_cropped <- asym_pnorm_cropped(max)
  #df <- cbind(df, data.frame(w_cropped))
  #df_clean <- df[df$w_cropped != 0, ]
  df$w_sim <- (2 - abs(df$wx - df$wy)) / 2
  df <- df[df$wx != 0 | df$wy != 0, ]
  x <- df$x
  y <- df$y
  w_sim <- df$w_sim
  #if (length(df_clean$x) > 1 & length(df$x) > 3) {
  if (length(df$x) > 5) {
    soft_cos_sim <-   sum(w_sim * x * y) / ( sqrt( sum(w_sim * x^2) ) * sqrt( sum(w_sim * y^2) ) )
  } else {
    soft_cos_sim <- NaN
  }
  return(soft_cos_sim)
}

bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')
##### bad strains are library strains to be removed (according to Amelie's file, and Hannes's instruction, those should be removed)

### gsp1 point mutant E-MAP data merged with the scaled SGA data
load("basic_E-MAP_data/SGA_Gsp1emap_merged.RData") 
ubermap <- ubermap %>% filter(! library_ORF %in% bad_strains)

library_strains_out_of_zero_box_for_wt_gsp1 <- ubermap %>% 
  filter(query_uniq == "GSP1-NAT" & (score < -3 | score > 2)) %>%
  unique()
##### 5 genes have significant (all POSITIVE!) scores with the WT Gsp1 (SAC3, PMP3, THP3, NUP133 and NUP188)
### all of them except for PMP3 are clearly relevant for nuclear transport so don't remove them from the ubermap
### however, when counting if there are enough library genes with significant scores to calculate a correlation, DON'T count those
  

## mutants only
ubermap.Gsp1_mutants.only <- ubermap %>% filter( query_ORF == "YLR293C")
ubermap.Gsp1_mutants_emap.only <- ubermap.Gsp1_mutants.only %>% filter(! grepl(pattern = "tsa", query_uniq))
all_gsp1_mutants <- ubermap.Gsp1_mutants.only %>% pull(query_uniq) %>% unique()
# save the mutant only data
write_tsv(ubermap.Gsp1_mutants.only, path = "basic_E-MAP_data/20181230_Gsp1mut_SGA_and_emap.txt")

ubermap.Gsp1_mutants.only %>% 
  ggplot(., mapping = aes(x = score)) + 
  geom_density() + ggtitle("Distribution of E-MAP scores in Gsp1 mutants")

#### look at the distribution of scores for Gsp1 mutants
ubermap.Gsp1_mutants.only %>% 
  ggplot(ubermap.mutants.only, 
         mapping = aes(x = score, 
             y = reorder(query_uniq, score, FUN = function(x) quantile(x, prob = 0.05, na.rm = T)))) +
         geom_boxplot() +
        ylab("Mutants ordered by the 5th E-MAP score percentile") +
        xlab("E-MAP score")

lib_genes_with_SGA_Gsp1 <- ubermap.Gsp1_mutants.only %>% 
  filter(query_uniq == "YLR293C_tsa256") %>% 
  pull(library_ORF) %>% unique()
ubermap.Gsp1_mutants.only_SGA_mut_lib <-  ubermap.Gsp1_mutants.only %>% 
  filter( library_ORF %in% lib_genes_with_SGA_Gsp1)

ubermap.per.mutant.95conf <- ubermap.Gsp1_mutants.only %>%
  group_by(query_uniq) %>%
  summarize(lower_quant = quantile(score, probs = .025, na.rm = T), 
            upper_quant = quantile(score, probs = .975, na.rm = T),
            upper_whisker = boxplot(score, plot = F)$stats[5, ],
            lower_whisker = boxplot(score, plot = F)$stats[1, ]) %>%
  mutate(category = factor(query_uniq, levels = query_uniq)) %>%
  arrange(lower_quant)
ubermap.per.mutant.95conf_strong <- ubermap.per.mutant.95conf %>% 
  filter(lower_quant < -3)
ubermap.per.mutant.95conf %>% 
  ggplot(., mapping = aes(x = lower_quant, y = upper_quant)) + geom_point() +
  geom_label_repel(data = ubermap.per.mutant.95conf_strong, aes(label = query_uniq)) +
  xlab("E-MAP score - lower 2.5 %") + ylab("E-MAP score - higher 97.5 %")
  
ubermap.per.mutant.95conf %>% 
  ggplot(., mapping = aes(x = lower_whisker, y = upper_whisker)) + geom_point() +
  geom_label_repel(data = ubermap.per.mutant.95conf_strong, aes(label = query_uniq)) +
  xlab("E-MAP score - lower whisker") + ylab("E-MAP score - upper whisker")

#### array genes per query
ubermap %>% 
  group_by(query_uniq) %>% 
  summarise("n_lib_genes" = n()) %>% 
  ggplot(aes(x = n_lib_genes)) + geom_histogram()

### queries per library gene
ubermap %>% 
  group_by(library_ORF) %>% 
  summarise("n_queries" = n()) %>% 
  filter(n_queries > 500) %>%
  ggplot(aes(x = n_queries)) + geom_histogram()
library_ORFs_to_keep <- ubermap %>% 
  group_by(library_ORF) %>% 
  summarise("n_queries" = n()) %>% 
  filter(n_queries > 500) %>% 
  pull(library_ORF) %>% unique()

#### array genes per query
ubermap %>% 
  filter(library_ORF %in% library_ORFs_to_keep) %>% 
  group_by(query_uniq) %>% 
  summarise("n_lib_genes" = n()) %>% 
  ggplot(aes(x = n_lib_genes)) + geom_histogram()
queries_to_keep <- ubermap %>% 
  filter(library_ORF %in% library_ORFs_to_keep) %>% 
  group_by(query_uniq) %>% 
  summarise("n_lib_genes" = n()) %>% 
  filter(n_lib_genes > 100) %>% 
  pull(query_uniq) %>% unique()


### skip this on 20190129
# ubermap <- ubermap %>%
#   filter(query_uniq %in% queries_to_keep & library_ORF %in% library_ORFs_to_keep)

#### pairwise plots of all Gsp1 mutants and core partners
core_partner_gene_names <- c("SRM1", "MOG1", "YRB1","YRB2", "YRB30", "RNA1", "LOS1", "MSN5", "MTR10",
                             "CSE1", "CRM1", "CEX1", "SRP1", "KAP95", "KAP120", "KAP123", "KAP104", 
                           "KAP114", "PSE1", "NMD5", "NUP42", "NTF2")
#core_partner_gene_names <- c("SRM1", "MOG1", "YRB1", "RNA1",
 #                           "CSE1", "CRM1", "KAP95", "PSE1", "NTF2")
core_ubermap <- ubermap %>% 
  filter(query_gene_name %in% core_partner_gene_names) %>% 
  separate(query_uniq, into = c("query_ORF", "pert"), sep = "_") %>% 
  mutate("partner" = str_c(query_gene_name, pert, sep = "_")) %>% 
  select(partner, library_ORF, "partner_score" = score)
Gsp1_mutants <- ubermap %>% filter(query_ORF == "YLR293C") %>% pull(query_uniq) %>% unique()
plot_lims <- range(ubermap.Gsp1_mutants.only$score, na.rm = T)
plots <- list()
correlations_with_partners <- tibble()

for (m in seq_along(Gsp1_mutants)) {
  mut <- Gsp1_mutants[m]
  mut_emap <- ubermap  %>% 
    filter(query_uniq == mut) %>% 
    select("Gsp1_mutant" = query_uniq, library_ORF, "Gsp1_score" = score) %>% 
    inner_join(., core_ubermap, by = "library_ORF") 
  correlations <- mut_emap %>%  
    group_by(partner) %>% 
    summarise("pearson" = cor(Gsp1_score, partner_score, use = "pairwise.complete.obs"),
              "w_pearson" = w_pearson(Gsp1_score, partner_score),
              "CosSim" = SoftCosSim_cropped(Gsp1_score, partner_score)
              ) %>% 
    mutate("subtitle" = str_c(partner, "cor =", round(w_pearson, 2), "CosSim =", round(CosSim, 2), sep = " ")) %>% 
    ungroup() %>% 
    inner_join(., mut_emap, by = "partner") %>% 
    arrange(desc(w_pearson))
  correlations_with_partners <- correlations %>% 
    select(partner, w_pearson, CosSim, Gsp1_mutant) %>% 
    unique() %>% 
    bind_rows(., correlations_with_partners)
  partners_ordered_by_corr <- correlations %>% pull(subtitle) %>% unique()
  plots[[m]] <- correlations %>% 
    mutate("subtitle" = factor(subtitle, partners_ordered_by_corr)) %>% 
    arrange(subtitle) %>% 
    ggplot(aes(x = Gsp1_score, y = partner_score)) + 
    geom_point() + facet_wrap(~ subtitle) + xlim(plot_lims) + ylim(plot_lims) +
    ggtitle(mut)
    
}
pdf("E-MAP_SGA_individual_scatterplots.pdf", height = 10, width = 18)
print(plots)
dev.off()






core_ubermap <- ubermap %>% 
  filter(query_gene_name %in% core_partner_gene_names) %>% 
  filter(library_ORF %in% lib_genes_with_SGA_Gsp1) %>% 
  separate(query_uniq, into = c("query_ORF", "pert"), sep = "_") %>% 
  mutate("partner" = str_c(query_gene_name, pert, sep = "_")) %>% 
  select(partner, library_ORF, "partner_score" = score)
plot_lims <- range(ubermap.Gsp1_mutants.only_SGA_mut_lib$score, na.rm = T)
plots <- list()
correlations_with_partners <- tibble()

for (m in seq_along(Gsp1_mutants)) {
  mut <- Gsp1_mutants[m]
  mut_emap <- ubermap.Gsp1_mutants.only_SGA_mut_lib  %>% 
    filter(query_uniq == mut) %>% 
    select("Gsp1_mutant" = query_uniq, library_ORF, "Gsp1_score" = score) %>% 
    inner_join(., core_ubermap, by = "library_ORF") 
  correlations <- mut_emap %>%  
    group_by(partner) %>% 
    summarise("pearson" = cor(Gsp1_score, partner_score, use = "pairwise.complete.obs"),
              "w_pearson" = w_pearson(Gsp1_score, partner_score),
              "CosSim" = SoftCosSim_cropped(Gsp1_score, partner_score)
    ) %>% 
    mutate("subtitle" = str_c(partner, "cor =", round(w_pearson, 2), "CosSim =", round(CosSim, 2), sep = " ")) %>% 
    ungroup() %>% 
    inner_join(., mut_emap, by = "partner") %>% 
    arrange(desc(w_pearson))
  correlations_with_partners <- correlations %>% 
    select(partner, w_pearson, CosSim, Gsp1_mutant) %>% 
    unique() %>% 
    bind_rows(., correlations_with_partners)
  partners_ordered_by_corr <- correlations %>% pull(subtitle) %>% unique()
  plots[[m]] <- correlations %>% 
    mutate("subtitle" = factor(subtitle, partners_ordered_by_corr)) %>% 
    arrange(subtitle) %>% 
    ggplot(aes(x = Gsp1_score, y = partner_score)) + 
    geom_point() + facet_wrap(~ subtitle) + xlim(plot_lims) + ylim(plot_lims) +
    ggtitle(mut)
  
}
pdf("E-MAP_SGA_individual_scatterplots_SGA_lib_only.pdf", height = 10, width = 18)
print(plots)
dev.off()



correlations_with_partners %>% 
  ggplot(aes(x = w_pearson)) + geom_histogram() + facet_wrap(~ partner) +
  ggtitle("Distribution of weighted Pearson correlations with Gsp1 mutants per regulator gene mutant")

correlations_with_partners %>% 
  ggplot(aes(x = CosSim)) + geom_histogram() + facet_wrap(~ partner) +
  ggtitle("Distribution of soft cosine similarities with Gsp1 mutants per regulator gene mutant")

selected_mutants_to_look_at <- c("T34E", "T34A", "T34Q", "T34G", "R108L", "R108A", "R108G", "R108Q", "R108I", 
                                "D79S", "D79A", "K143W", "H141E", "R78K", "Y157A", "K101R", "YLR293C_tsa256", "Q147E",
                                "R112S", "G80A")

correlations_with_partners %>% 
  filter(Gsp1_mutant %in% selected_mutants_to_look_at) %>% 
  ggplot(aes(x = partner, y = Gsp1_mutant, size = w_pearson, color = w_pearson)) + geom_point()  +
  scale_color_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Weighted pearson corr with regulator gene mutants for selected Gsp1 mutants")


sga <- ubermap %>% 
  mutate("query_uniq" = ifelse(query_uniq == "NTER3XFLAG_WT", "NTER3XFLAG-WT", query_uniq)) %>% 
  mutate("query_uniq" = ifelse(query_uniq == "CTER3XFLAG_WT", "CTER3XFLAG-WT", query_uniq)) %>% 
  filter(grepl(query_uniq, pattern = "_")) %>% 
  separate(query_uniq, into = c("query_ORF", "pert"), extra = "merge", sep = "_") %>% 
  mutate("query_uniq" = str_c(query_gene_name, pert, sep = "_")) %>% 
  select(-pert)
Gsp1_mutants <- ubermap %>% 
  filter(! grepl(query_uniq, pattern = "_"))

ubermap <- bind_rows(Gsp1_mutants, sga)

save(ubermap, file = "basic_E-MAP_data/20190129_ubermap_for_correlations.RData")






corr_with_partners_spread <- correlations_with_partners %>% 
  select(partner, Gsp1_mutant, CosSim) %>% 
  spread(key = partner, value = CosSim)
corr_with_partners_spread <- data.frame(corr_with_partners_spread)
dat <- corr_with_partners_spread[, 2:30]  # numerical columns
rownames(dat) <- corr_with_partners_spread$Gsp1_mutant
row.order <- hclust(dist(dat))$order # clustering
col.order <- hclust(dist(t(dat)))$order
dat_new <- dat[row.order, col.order] # re-order matrix accoring to clustering

df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
names(df_molten_dat) <- c("Gsp1_mutant", "Partner_protein", "correlation")

ggplot(data = df_molten_dat,
       aes(y = Partner_protein, x = Gsp1_mutant, color = correlation, size = correlation)) + 
  geom_point() +
  scale_color_gradient2() +
  #scale_fill_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Weighted Pearson correlation between Gsp1 mutants and TS mutants pf partners from SGA data")

