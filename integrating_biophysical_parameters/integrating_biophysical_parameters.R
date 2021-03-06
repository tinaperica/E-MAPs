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
GEF_kin <- read_tsv("integrating_biophysical_parameters/GEF_kinetics_param.txt", col_names = T) %>% 
  gather("measure", "value", -mutant)
GAP_kin <- read_tsv("integrating_biophysical_parameters/GAP_MM_parameters_clean.txt", col_names = T) %>% 
  gather("measure", "value", -mutant)
NMR_data <- read_tsv("integrating_biophysical_parameters/20190224_nmr_data.txt", col_names = T) %>% 
  select("mutant" = "gsp1_mutation", "gamma2" = gamma2_percent) %>% 
  gather("measure", "value", -mutant)
apms_foldchange <- read_tsv("integrating_biophysical_parameters/apms_log2_fold_change.txt", col_names = T) %>% 
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
### input APMS based GAP GEF ratio
apms_GAP_GEF_diff <- read_tsv("integrating_biophysical_parameters/tag_averagged_gap_minus_gef_ln_APMS_fold_change_from_WT.txt")

### hierarchical clustering of the emap data
# spread_emap <- e.map %>% 
#   filter(mutant %in% unique(all_parameters$mutant)) %>% 
#   select(library_gene_name, mutant, score) %>% 
#   spread(library_gene_name, score)
# spread_emap <- data.frame(spread_emap)
# rownames(spread_emap) <- spread_emap$mutant
# spread_emap <- spread_emap[, -1]
# spread_emap <- Filter(function (x) !all (is.na(x)), spread_emap)  # this removes columns that are all NA
# head(spread_emap[1:10, 1:10])
# hclust_mut_data <- hcut(spread_emap, k = 4)
# order_mutants <- hclust_mut_data$labels[hclust_mut_data$order]
# fviz_dend(hclust_mut_data, k_colors = brewer.pal(4, "Set1"),
#           main = "mutants clustered by E-MAP score, k = 4")  # from ‘factoextra’ package 
# 
# 
# 
# 
# cluster_by_parameters <- function(emap, biophy_params, k) {
#   plots <- list()
#   param_subset <- all_parameters %>% 
#     filter(measure %in% biophy_params)
#   mutants_with_all_params <- param_subset %>% 
#     select(mutant, measure) %>% 
#     group_by(mutant) %>% 
#     summarise("n_param" = n()) %>% 
#     ungroup() %>% 
#     mutate("max_n_param" = max(n_param)) %>% 
#     filter(n_param == max_n_param) %>% 
#     pull(mutant)
#   emap.subset <- emap %>% 
#     filter(mutant %in% mutants_with_all_params)
#   
#   ### hierarchical clustering of the emap data
#   spread_emap <- emap.subset %>% 
#     #filter(mutant %in% unique(all_parameters$mutant)) %>% 
#     select(library_gene_name, mutant, score) %>% 
#     spread(library_gene_name, score)
#   spread_emap <- data.frame(spread_emap)
#   rownames(spread_emap) <- spread_emap$mutant
#   spread_emap <- spread_emap[, -1]
#   spread_emap <- Filter(function (x) !all (is.na(x)), spread_emap)  # this removes columns that are all NA
#   head(spread_emap[1:10, 1:10])
#   hclust_emap_data <- hcut(spread_emap, k = k)
#   
#   mut_groups <- data.frame(cutree(hclust_emap_data, k = k))
#   names(mut_groups) <- "group"
#   mut_groups <- tibble("mutant" = rownames(mut_groups), "group" = mut_groups$group) %>% 
#     inner_join(., tibble("group" = seq(1, k, 1), "color" = brewer.pal(k, "Set1")), by = "group") %>% 
#     arrange(group)
#   order_mutants <- hclust_emap_data$labels[hclust_emap_data$order]
#   mut_groups <- mut_groups %>% 
#     mutate("mut_fact" = factor(mutant, order_mutants)) %>% 
#     arrange(mut_fact)
#   plots[["by_emap"]] <- fviz_dend(hclust_emap_data, label_cols =  mut_groups$color, k_colors = "black",
#             main = str_c("mutants clustered by E-MAP score, k = ", k))  # from ‘factoextra’ package 
#   
#   #### cluster by parameters
#   param_subset_spread <- param_subset %>% 
#     filter(mutant %in% mutants_with_all_params) %>% 
#     select(-raw_value) %>% 
#     spread(measure, value)
#   param_subset_spread <- data.frame(param_subset_spread)
#   rownames(param_subset_spread) <- param_subset_spread$mutant
#   param_subset_spread <- param_subset_spread[, -1]
#   hclust_param_data <- hcut(param_subset_spread, k = k)
#   order_mutants_by_param <- hclust_param_data$labels[hclust_param_data$order]
#   mut_groups <- mut_groups %>% 
#     mutate("mut_fact" = factor(mutant, order_mutants_by_param)) %>% 
#     arrange(mut_fact)
#   plots[["by_param"]] <- fviz_dend(hclust_param_data, label_cols =  mut_groups$color, k_colors = "black", 
#             main = str_c( c("mutants clustered by: (", biophy_params, ") -> k =", k), collapse = " "))
#   pdf(file = str_c(c("integrating_biophysical_parameters/", date, "_", biophy_params, ".pdf"), collapse = ""), width = 10)
#   print(plots)
#   dev.off()
# }
# 
# available_parameters <- all_parameters %>% pull(measure) %>% unique()
# cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km"), k = 6)
# cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat_Km", "GEF_kcat_Km"), k = 4)
# 
# cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2"), k = 4)
# 
# cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1"), k = 4)
# cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "MOG1"), k = 4)
# cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "YRB1"), k = 4)
# cluster_by_parameters(emap = e.map, biophy_params = c("SRM1", "RNA1", "MOG1", "YRB1"), k = 4)
# cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2", "SRM1", "RNA1"), k = 4)
# cluster_by_parameters(emap = e.map, biophy_params = c("GAP_kcat", "GAP_Km", "GEF_kcat", "GEF_Km", "gamma2", "SRM1", "RNA1", "MOG1", "YRB1"), k = 4)



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
  ggplot(aes(x = gamma2, y = log(GAP_kcat_Km))) +
  geom_point(size = 5, alpha = 0.5)  +
  geom_text_repel(aes(label = mutant), size = 7) +
  ylab("GAP hydrolysis relative kcat/Km\n\nln(kcat/Km(MUT) / kcat/Km(WT))\n") +
  xlab("\n% in hydrolysis ready conformation") +
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18)
        )
ggsave(filename = "integrating_biophysical_parameters/GAP_enzyme_eff_vs_gamma2.pdf", width = 9, height = 5)


relative_parameters <- all_parameters %>% 
  filter(measure %in% c("GEF_kcat", "GEF_Km", "GEF_kcat_Km", 
                        "GAP_kcat", "GAP_Km", "GAP_kcat_Km", "gamma2")) %>% 
  select(-value) %>% 
  inner_join(., WT_parameters, by = c("measure")) %>% 
  mutate("rel_value" = raw_value/wt_value) %>% 
  select(-raw_value, -wt_value) 

rel_GAP_GEF_efficiency <- all_parameters %>% 
  filter(measure %in% c("GEF_kcat_Km", "GAP_kcat_Km")) %>% 
  select(-value) %>% 
  inner_join(., WT_parameters, by = c("measure")) %>% 
  mutate("rel_value" = raw_value/wt_value) %>% 
  select(-raw_value, -wt_value) 
  

rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  ggplot(aes(x = log(GEF_kcat_Km), y = log(GAP_kcat_Km), label = mutant)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
  geom_point() +
  geom_text_repel() +
  xlab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GEF mediated nucleotide exchange") +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  xlim(c(-6, 1)) + ylim(c(-6, 1)) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

ggsave(filename = "integrating_biophysical_parameters/GAP_enzyme_eff_vs_GEF_enzyme_eff.pdf", width = 9, height = 7)

rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  inner_join(., apms_GAP_GEF_diff) %>% 
  ggplot(aes(x = log(GEF_kcat_Km), y = log(GAP_kcat_Km), label = mutant)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
  geom_point(aes(color = gap_minus_gef_FC), size = 5) +
  geom_text_repel() +
  scale_color_gradient2() +
  labs(color = "GAP - GEF\nln(fold change MUT/WT)") +
  xlab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GEF mediated nucleotide exchange") +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  xlim(c(-6, 1)) + ylim(c(-6, 1)) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("integrating_biophysical_parameters/GAP_enzyme_eff_vs_GEF_enzyme_eff_color_by_APMS_GAP_GEF.pdf", width = 10, height = 7)

# same plot, but just GAP APMS log2FC, then just GEF APMS log2FC
apms_GAP_GEF <- read_tsv("integrating_biophysical_parameters/tag_averaged_apms_GAP_GEF_log2FC.txt")

rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  inner_join(., apms_GAP_GEF) %>% 
  select(-SRM1) %>% 
  ggplot(aes(x = log(GEF_kcat_Km), y = log(GAP_kcat_Km), label = mutant)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
  geom_point(aes(color = RNA1), size = 5) +
  geom_text_repel() +
  scale_color_gradient2() +
  labs(color = "GAP ln(fold change MUT/WT)") +
  xlab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GEF mediated nucleotide exchange") +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  xlim(c(-6, 1)) + ylim(c(-6, 1)) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("integrating_biophysical_parameters/GAP_enzyme_eff_vs_GEF_enzyme_eff_color_by_APMS_GAP_only.pdf", width = 10, height = 7)

rel_GAP_GEF_efficiency %>% 
  spread(measure, rel_value) %>% 
  inner_join(., apms_GAP_GEF) %>% 
  select(-RNA1) %>% 
  ggplot(aes(x = log(GEF_kcat_Km), y = log(GAP_kcat_Km), label = mutant)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
  geom_point(aes(color = SRM1), size = 5) +
  geom_text_repel() +
  scale_color_gradient2() +
  labs(color = "GEF ln(fold change MUT/WT)") +
  xlab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GEF mediated nucleotide exchange") +
  ylab("ln(kcat/Km(MUT) / kcat/Km(WT)) of GAP mediated GTP hydrolysis") +
  xlim(c(-6, 1)) + ylim(c(-6, 1)) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("integrating_biophysical_parameters/GAP_enzyme_eff_vs_GEF_enzyme_eff_color_by_APMS_GEF_only.pdf", width = 10, height = 7)


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
  ylab("ln(relative GAP / relative GEF efficiency)") +
  xlab("Gsp1 point mutant") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
        axis.title = element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", alpha = 0.5)

ggsave(filename = "integrating_biophysical_parameters/titration_of_relative_GAP_and_GEF_efficiency.pdf", width = 10)
write_tsv(GAP_GEF_ratio, "titration_curves/GAP_GEF_ratio.txt")

mut_ranked_by_GAP_eff <- GAP_GEF_ratio %>% arrange(GAP_kcat_Km) %>% 
  select(mutant, GAP_kcat_Km)
write_tsv(mut_ranked_by_GAP_eff, "integrating_biophysical_parameters/ranked_by_GAP_eff.txt")

mut_ranked_by_GEF_eff <- GAP_GEF_ratio %>% arrange(GEF_kcat_Km) %>% 
  select(mutant, GEF_kcat_Km)
write_tsv(mut_ranked_by_GEF_eff, "integrating_biophysical_parameters/ranked_by_GEF_eff.txt")



GAP_GEF_ratio <- relative_parameters %>% 
  spread(measure, rel_value) %>% 
  mutate("GAP/GEF efficiency" = GAP_kcat_Km / GEF_kcat_Km,
         "GAP kcat / GEF kcat" = GAP_kcat / GEF_kcat) %>% 
  filter(! is.na(`GAP/GEF efficiency`)) %>% 
  arrange(`GAP/GEF efficiency`) 
mutants_by_GAP_GEF_ratio <- GAP_GEF_ratio %>% 
  pull(mutant) %>% unique()

GAP_GEF_ratio %>% 
  #select(mutant, `GAP/GEF efficiency`, `GAP kcat / GEF kcat`) %>% 
  gather(parameter, value, -mutant) %>% 
  mutate("mutant" = factor(mutant, mutants_by_GAP_GEF_ratio)) %>% 
  ggplot(aes(x = mutant, y = log(value), color = parameter, group = 1)) +
  geom_point(size = 6, alpha = 0.5) +
  ylab("ln(value)") +
  xlab("Gsp1 point mutant") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
        axis.title = element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", alpha = 0.5)

write_tsv(GAP_GEF_ratio, "titration_curves/GAP_GEF_ratio_all.txt")

GAP_GEF_ratio %>% 
  select(mutant, gamma2, GAP_kcat_Km) %>% 
  gather(parameter, value, -mutant) %>% 
  mutate("mutant" = factor(mutant, mutants_by_GAP_GEF_ratio)) %>% 
  ggplot(aes(x = mutant, y = log(value), color = parameter, group = 1)) +
  geom_point(size = 6, alpha = 0.5) +
  ylab("ln(value)") +
  xlab("Gsp1 point mutant") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
        axis.title = element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", alpha = 0.5)
