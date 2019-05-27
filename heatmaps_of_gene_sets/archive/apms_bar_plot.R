# Code to make bar plot for titration curves with bars

# Read in APMS data for tag barplot
GAP_GEF_log2FC <- read_tsv("datasets/apms_log2_fold_change.txt") %>%
  filter(gene_name %in% c("SRM1", "RNA1") & norm == "eqM")

tag_avrg <- GAP_GEF_log2FC %>%
  arrange(mutant, tag, gene_name) %>%
  group_by(mutant, gene_name) %>%
  summarise("tag_avg_log2FC" = mean(log2FC, na.rm = T),
            "tag_avg_log2FC_sd" = sd(log2FC, na.rm = T)) %>%
  ungroup() %>%
  mutate("tag_avg_log2FC_sd" = ifelse(is.na(tag_avg_log2FC_sd), 0, tag_avg_log2FC_sd)) %>%
  inner_join(., GAP_GEF_log2FC, by = c("mutant", "gene_name")) %>%
  select(mutant, gene_name, tag_avg_log2FC, tag_avg_log2FC_sd, tag)

mut_with_both_tags <- tag_avrg %>%
  group_by(mutant, gene_name) %>%
  summarize("count" = n()) %>%
  filter(count > 1) %>% pull(mutant) %>% unique()

tag_avrg <- tag_avrg %>%
  mutate("tag"= ifelse(is.element(mutant, mut_with_both_tags), "both", tag)) %>% unique()

tag_avrg_gap_minus_gef <- tag_avrg %>%
  select(mutant, tag, everything()) %>%
  gather(variable, value, -(mutant:gene_name)) %>%
  unite(temp, gene_name, variable) %>%
  spread(temp, value) %>%
  group_by(mutant, tag) %>%
  summarize("gap_minus_gef_FC" = RNA1_tag_avg_log2FC - SRM1_tag_avg_log2FC,
            "gap_minus_gef_sd" = sqrt( (RNA1_tag_avg_log2FC_sd)^2 + (SRM1_tag_avg_log2FC_sd)^2)) %>%
  ungroup()

apms_with_NA_bar <- left_join(unique(select(emap, mutant)), tag_avrg_gap_minus_gef, by = 'mutant')


# AP-MS bar plot
apms_plot_bar <-
  apms_with_NA_bar %>%
  mutate("mutant" = factor(mutant, ordered_genes)) %>%
  ggplot(aes(x = mutant, y = gap_minus_gef_FC, fill = tag)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin = gap_minus_gef_FC - gap_minus_gef_sd,
                    ymax = gap_minus_gef_FC + gap_minus_gef_sd), width = 0.5) +
  theme_grey() +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9),
        legend.direction = 'horizontal',
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_rect(size = 0.1)) +
  ylab("GAP - GEF\ndelta log2FC") + xlab(' ')