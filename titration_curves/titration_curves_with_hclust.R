#### this script makes "titration curves" of E-MAP scores for Gsp1 mutants
#### curves are made per Wodak complex
### load in E-MAP data
library(tidyverse)
abs_max <- function(x) {
  x <- x[! is.na(x)]
  if (length(x) > 1) {
    abs_max <- x[which.max( abs(x))] 
  } else if (length(x) == 1) {
    abs_max <- x
  } else {
    abs_max <- NaN
  }
  return(abs_max)
}
load("basic_E-MAP_data/20180630_emap_data_for_corr_all_clusters.RData")
genes <- c("GSP1", "MOG1", "RNA1", "SRM1", "YRB1", "PSE1")
emap.data <- as_tibble(ubermap$ubermap) %>% 
  filter(Gene.gene_name %in% genes) %>% 
  filter(Gene_uniq != "CTER3XFLAG WT" & Gene_uniq != "NTER3XFLAG WT") %>% 
  ungroup() %>% 
  mutate(Gene_uniq = ifelse(Gene.gene_name != "GSP1", Gene.gene_name, Gene_uniq)) %>% 
  select(-Gene.gene_name, -ORF)
library_clusters <- as_tibble(ubermap$library_clusters) %>% 
  rename("library.ORF" = ORF, "library.gene_name" = gene_name)
clusters <- ubermap$clusters
rm(ubermap)
wodak_complexes <- read_tsv("titration_curves/all_wodak_complexes.txt", col_names = T) %>% 
  select("library.ORF" = ORF, "library.gene_name" = Name, "cluster" = Complex)
cmplxs <- wodak_complexes %>% pull(cluster) %>% unique()
clusters_and_complexes <- bind_rows(wodak_complexes, library_clusters)
clst_cmplxs <- clusters_and_complexes %>% pull(cluster) %>% unique()
merged_data <- emap.data %>% 
  inner_join(., clusters_and_complexes) 

choose_to_plot <- merged_data %>% 
  select(Gene_uniq, library.ORF, library.gene_name, score, random_high_score, cluster) %>% 
  filter(! is.na(score)) %>% 
  group_by(Gene_uniq, cluster) %>% 
  summarise("stdev_score" = sd(score, na.rm = T), 
            "stdev_random_high_score" = sd(random_high_score, na.rm = T),
            "max_abs_score" = abs_max(score), 
            "max_abs_random_high_score" = abs_max(random_high_score),
            "abs_max_score" = abs(max(score, na.rm = T)),
            "abs_max_random_high_score" = abs(max(random_high_score, na.rm = T))) %>% 
  filter(! is.na(abs_max_score))
clusters_to_plot <- choose_to_plot %>% 
  #filter(stdev_random_high_score > stdev_score & abs_max_score > abs_max_random_high_score) %>% 
  #filter(stdev_random_high_score > stdev_score) %>% 
  pull(cluster) %>% unique()
clusters_to_plot <- c("COMA complex", "Elongator complex")
#### print titration curves by wodak complex/cluster chosen for plotting
y_lim <- c(min(c(emap.data$score), na.rm = T), max(c(emap.data$score), na.rm = T))
plots <- list()
plot_counter <- 1
hclust_plots <- list()
hclust_plot_counter <- 1
heatmap_plots <- list()
heatmap_counter <- 1
for (i in seq_along(clusters_to_plot)) {
  clust <- clusters_to_plot[i]
  all_cluster_complex_genes <- clusters_and_complexes %>% 
    filter(cluster == clust) %>% pull(library.gene_name) %>% unique()
  temp <- merged_data %>% 
    filter(cluster == clust) %>% 
    select(Gene_uniq, score, "library gene" = library.gene_name, cluster)
  cluster_complex_genes_with_data <- temp %>% 
    pull(`library gene`) %>% unique()
  if (length(cluster_complex_genes_with_data) > 1) {
    missing_data_genes <- all_cluster_complex_genes[! all_cluster_complex_genes %in% cluster_complex_genes_with_data] 
    temp <- temp %>% bind_rows(., tibble("Gene_uniq" = "A180T", "score" = NA, 
                                       "library gene" = missing_data_genes, "cluster" = clust))
    temp <- temp %>% complete(Gene_uniq, `library gene`, cluster)
    temp <- temp %>% 
      group_by(`library gene`) %>% 
      mutate("abs_max" = abs(abs_max(score))) %>% 
      ungroup()
    mutants_spread <- temp %>% 
      select(Gene_uniq, `library gene`, score) %>% 
      mutate("score" = ifelse(is.na(score), 0, score)) %>% 
      spread(`library gene`, score)
    mutants_ordered_n <- hclust(dist(mutants_spread))$order
    mutants_ordered_by_hclust <- temp %>%
      pull(Gene_uniq) %>% unique()
    mutants_ordered_by_hclust <- mutants_ordered_by_hclust[mutants_ordered_n]
    lib_genes_spread <- temp %>% 
      select(Gene_uniq, `library gene`, score) %>% 
      mutate("score" = ifelse(is.na(score), 0, score)) %>% 
      spread(Gene_uniq, score )
    lib_genes_ordered_n <- hclust(dist(lib_genes_spread))$order
    lib_genes_ordered_by_hclust <- temp %>%
      pull(`library gene`) %>% unique()
    lib_genes_ordered_by_hclust <- lib_genes_ordered_by_hclust[lib_genes_ordered_n]
    mutants_ordered_by_mean_score <- temp %>%
       group_by(Gene_uniq) %>% 
       summarise("mean_score" = mean(score, na.rm = T)) %>% 
       arrange(mean_score) %>% pull(Gene_uniq)
    temp <- temp %>%
      mutate("mutant" = factor(Gene_uniq, mutants_ordered_by_mean_score)) %>% 
      arrange(mutant)
    hclust_temp <- temp %>% 
      mutate("mutant" = factor(Gene_uniq, mutants_ordered_by_hclust),
             "library gene" = factor(`library gene`, lib_genes_ordered_by_hclust)) %>% 
      arrange(mutant, `library gene`)
    if (length(cluster_complex_genes_with_data) < 10) {
      # hclust_plots[[hclust_plot_counter]] <- hclust_temp %>%  ggplot(., aes(x = mutant, y = score)) + 
      #   geom_abline(intercept = 0, slope = 0, col = "gray") +
      #   geom_point(aes(color = `library gene`), size = 3) + 
      #   theme_bw() +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   ggtitle(clust) + ylab("E-MAP score") +
      #   ylim(y_lim[1], y_lim[2])
      # hclust_plot_counter <- hclust_plot_counter + 1
      plots[[plot_counter]] <- temp %>%  ggplot(., aes(x = mutant, y = score)) + 
        geom_abline(intercept = 0, slope = 0, col = "gray") +
        geom_point(aes(color = `library gene`), size = 3) + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(clust) + ylab("E-MAP score") +
        ylim(y_lim[1], y_lim[2])
      plot_counter <- plot_counter + 1
    }
    # if (length(cluster_complex_genes_with_data) > 5) {
    #   heatmap_plots[[heatmap_counter]] <- hclust_temp %>% 
    #     ggplot(., aes(x = mutant, y = `library gene`)) +
    #     geom_tile(aes(fill = score)) +
    #      scale_fill_gradientn(colours = c("#E69F00", "#F0E442", "#F0E442", "black","black", "black", "#56B4E9", "#0072B2"),
    #                           values = rescale(c(y_lim[2], 4, 2.5, 1, 0, -2, -5, y_lim[1])),
    #                           guide = "colourbar") +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #     ggtitle(clust)
    #   heatmap_counter <- heatmap_counter + 1
    # }
  }
}

pdf("titration_curves/titration_curves_COMA_and_elongator.pdf", width = 12)
print(plots)
dev.off()

pdf("titration_curves/titration_curves_clusters_and_complexes.pdf", width = 12)
print(plots)
dev.off()

pdf("titration_curves/hclust_titration_curves_clusters_and_complexes.pdf", width = 12)
print(hclust_plots)
dev.off()
pdf("titration_curves/E-MAP_score_heatmaps_clusters_and_complexes.pdf", height = 22, width = 12)
print(heatmap_plots)
dev.off()

#### make titration curves with interface information
emap.mut.only <- merged_data %>% 
  filter(! Gene_uniq %in% c("GSP1-NAT", "NTER3XFLAG WT", "CTER3XFLAG WT", "MOG1", "RNA1", "SRM1", "YRB1", "PSE1")) %>% 
  mutate("residue" = as.numeric(substring(Gene_uniq, 2, (nchar(Gene_uniq)-1))))
residues <- emap.mut.only %>% 
  pull(residue) %>% unique() %>% sort()
contacts_and_interfaces <- read_tsv("~/Documents/Gsp1_bioinformatics/Ran_structures/4A_contacts_and_interfaces.txt")
contacts_and_interfaces <- contacts_and_interfaces %>%
  filter(yeastresnum %in% residues) %>% 
  select(partner, "residue" = yeastresnum, mean_n_contacts)
partner_lst_df <- tibble("res" = integer(), "partner_proteins" = character())
for (res in residues) {
  temp <- contacts_and_interfaces %>% 
    filter(residue == res)
  partner_proteins <- temp %>% 
    arrange(mean_n_contacts) %>% pull(partner) %>% 
    str_c(collapse = " ")
  partner_lst_df <- add_row(partner_lst_df, res, partner_proteins)
}
emap.mut_merge <- emap.mut.only %>% 
  inner_join(., partner_lst_df, by = c("residue" = "res")) %>% 
  mutate("mutant" = str_c(Gene_uniq, partner_proteins, sep = " - "))



hclust_plots <- list()
hclust_plot_counter <- 1
for (i in seq_along(clusters_to_plot)) {
  clust <- clusters_to_plot[i]
  all_cluster_complex_genes <- clusters_and_complexes %>% 
    filter(cluster == clust) %>% pull(library.gene_name) %>% unique()
  temp <- emap.mut_merge %>% 
    filter(cluster == clust) %>% 
    select(mutant, score, "library gene" = library.gene_name, cluster)
  cluster_complex_genes_with_data <- temp %>% 
    pull(`library gene`) %>% unique()
  if (length(cluster_complex_genes_with_data) > 1) {
    missing_data_genes <- all_cluster_complex_genes[! all_cluster_complex_genes %in% cluster_complex_genes_with_data] 
    temp <- temp %>% bind_rows(., tibble("mutant" = "A180T - KAP95 YRB2 YRB1", "score" = NA, 
                                         "library gene" = missing_data_genes, "cluster" = clust))
    temp <- temp %>% complete(mutant, `library gene`, cluster)
    temp <- temp %>% 
      group_by(`library gene`) %>% 
      mutate("abs_max" = abs(abs_max(score))) %>% 
      ungroup()
    mutants_spread <- temp %>% 
      select(mutant, `library gene`, score) %>% 
      mutate("score" = ifelse(is.na(score), 0, score)) %>% 
      spread(`library gene`, score)
    mutants_ordered_n <- hclust(dist(mutants_spread))$order
    mutants_ordered_by_hclust <- temp %>%
      pull(mutant) %>% unique()
    mutants_ordered_by_hclust <- mutants_ordered_by_hclust[mutants_ordered_n]
    lib_genes_spread <- temp %>% 
      select(mutant, `library gene`, score) %>% 
      mutate("score" = ifelse(is.na(score), 0, score)) %>% 
      spread(mutant, score )
    lib_genes_ordered_n <- hclust(dist(lib_genes_spread))$order
    lib_genes_ordered_by_hclust <- temp %>%
      pull(`library gene`) %>% unique()
    lib_genes_ordered_by_hclust <- lib_genes_ordered_by_hclust[lib_genes_ordered_n]
    hclust_temp <- temp %>% 
      mutate("mut" = factor(mutant, mutants_ordered_by_hclust),
             "library gene" = factor(`library gene`, lib_genes_ordered_by_hclust)) %>% 
      arrange(mut, `library gene`)
    if (length(cluster_complex_genes_with_data) < 10) {
      hclust_plots[[hclust_plot_counter]] <- hclust_temp %>%  ggplot(., aes(x = mut, y = score)) + 
        geom_abline(intercept = 0, slope = 0, col = "gray") +
        geom_point(aes(color = `library gene`), size = 3) + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(clust) + ylab("E-MAP score") +
        ylim(y_lim[1], y_lim[2])
      hclust_plot_counter <- hclust_plot_counter + 1
    }
  }
}

pdf("titration_curves/hclust_titration_curves_clusters_and_complexes_partner_proteins_listed.pdf", width = 12)
print(hclust_plots)
dev.off()



