library(tidyverse)
library(ggrepel)
library(ggcorrplot)
library(gmp)
library(GGally)
bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')
##### bad strains are library strains to be removed (according to Amelie's file, and Hannes's instruction, those should be removed)


### hannes_averaged4 E-MAP score data for Gsp1 mutants
#### read it in and make it tidy
e.map <- read_tsv("basic_E-MAP_data/avg_merged_June2016_screen_for_Gia.txt", col_names = T)   ## use export for Gia because it has ORF names for library
e.map <- e.map %>% gather(library, score, -Gene) %>%
  filter(! library %in% bad_strains) %>%
  arrange(desc(Gene), desc(library)) %>%
  mutate("Gene" = gsub(Gene, pattern = "GSP1:", replacement = "", perl = T))

#### these are all the gsp1 mutants screened (including WT-NAT and 3xFLAG constructs)
### 59 query constructs all together
(Gsp1_queries <- unique(e.map[["Gene"]]))  
### ORF to gene name annotation from SGD
orf_gene_name_index <- read_delim("orf_gene_GO_sgd_annotation.txt", col_names = F, delim = "\t")
orf_index <- unique(tibble("orf" = orf_gene_name_index$X1, "gene_name" = orf_gene_name_index$X2))
rm(orf_gene_name_index)
#############################

### Gsp1 mutants merged with ubermap (merged with ubermap if there is an at least 100 library genes overlap between 
# chromatin library and the library against which that query was screened)
# (gene names) - this is only necessary to distinguish point mutants (they are all the same in the ORF file)
gene_names.ubermap <- read_delim("basic_E-MAP_data/20180108_gene_names_merge_w_Ubermap_100.txt", 
                           delim = "\t", col_names = T, n_max = length(Gsp1_queries))
gene_names.ubermap <- gene_names.ubermap %>% gather(library.ORF, score, -Gene) %>%
  mutate("Gene" = gsub(Gene, pattern = "GSP1 - ", replacement = "", perl = T))
mutants <- unique(gene_names.ubermap$Gene)
rm(gene_names.ubermap)
### Gsp1 mutants merged with ubermap (orf names)
# replace Gsp1 ORFs with mutation names (because they are all the same orf!)
ubermap <- read_delim("basic_E-MAP_data/20180108_orf_names_merge_w_Ubermap_100.txt", col_names = T, delim = "\t")
ubermap[ 1:length(mutants), 1 ] <- mutants

#### can't gather before sorting out those few Genes that occur multiple times (add a unique identifier to them)
### first find those Genes that are not unique (probably data from separate experiments/labs for the same thing? 
## Stll better not to mix them up)
(not_uniq <- ubermap %>% 
    group_by(Gene) %>% 
    summarize("count" = n()) %>% 
    filter(count > 1) %>% 
    ungroup()
)

#### check correlations between non-unique samples
not_uniq_uber <- ubermap %>% 
  filter(Gene %in% not_uniq$Gene) %>% 
  arrange(Gene) %>% 
  mutate("Gene" = str_c(Gene, 1:24, sep = "_"))
not_uniq_matrix <- as.matrix(not_uniq_uber[, -1])
rownames(not_uniq_matrix) <- not_uniq_uber$Gene
cormat <- round(cor(t(not_uniq_matrix[, -1]), use = "pairwise.complete.obs"), 2)
####### correlations between duplicate data doesn't look that good
cormat %>% ggcorrplot(hc.order = TRUE, outline.col = "white",insig = "blank", tl.cex = 7) + 
  ggtitle("Pearson correlations between repeated screens") 
### these correlations look really bad
## check if this is due to the small overlap in libraries?

not_uniq_uber_gather <- not_uniq_uber %>% 
  gather(library.ORF, score, -Gene) %>% 
  mutate("Gene.ORF" = gsub(Gene, pattern = " - .*", replacement = "", perl = T))
non_uniq_gene_ORFS <- not_uniq_uber_gather %>% pull(Gene.ORF) %>% unique()
plots <- list()
for (i in seq_along(non_uniq_gene_ORFS)) {
  plots[[i]] <- not_uniq_uber_gather %>% 
    filter(Gene.ORF == non_uniq_gene_ORFS[i]) %>%
    select(-Gene.ORF) %>% 
    spread(Gene, score) %>%
    select(-library.ORF) %>% 
    ggpairs()
}
print(plots) 
##### add a number in front of all the non-unique Gene
for (temp.gene in not_uniq$Gene) {
  row.names <- which(ubermap$Gene == temp.gene)
  for (i in seq_along(row.names)) {
    ubermap$Gene[row.names[i]] <- str_c(i, ubermap$Gene[row.names[i]], sep = "_")
  }
}
### now those few non-unique Genes look like this: 1_YDL155W - YDL155W
#### check
filter(ubermap, Gene == "1_YDL155W - YDL155W" | Gene == "2_YDL155W - YDL155W")
(not_uniq <- ubermap %>% group_by(Gene) %>% 
    summarize("count" = n()) %>% filter(count > 1))


# now I can gather the data
ubermap <- ubermap %>% gather(library.ORF, score, -Gene)

### for RevCorofCorr I will only keep the query genes that have at least a 300 overlap with 
### the library Gsp1 was screened against
queries_to_keep <- ubermap %>% 
  filter(! is.na(score)) %>% 
  group_by(Gene) %>% 
  mutate("libORFcount" = n()) %>% 
  ungroup() %>% 
  arrange(desc(libORFcount)) %>%
  filter(libORFcount > 500) %>% 
  pull(Gene) %>% unique()     
### total number of queries is 4926, 
# after filtering for a 500 library genes overlap that number is 4612

ubermap <- ubermap %>% 
  filter(Gene %in% queries_to_keep)
## mutants only
ubermap.mutants.only <- ubermap %>% filter( Gene %in% mutants) %>%
  #mutate("library.ORF" = gsub(library.ORF, pattern = "\\.", replacement = "-", perl = T)) %>% 
  ### this is to replace the YOR298C.A with YOR298C-A (R puts . instead of - when loading library genes as column names)
  ### deprecated step now that we have tidyverse!
  add_column("ORF" = "YLR293C") %>% 
  rename("Gene_uniq" = Gene) %>% 
  mutate("Gene.gene_name" = "GSP1") %>% 
  inner_join(., orf_index, by = c("library.ORF" = "orf")) %>%
  rename("library.gene_name" = gene_name)

# save the mutant only data
write_tsv(ubermap.mutants.only, path = "basic_E-MAP_data/20180920_preprocessed_ubermap_mut_only.txt")

ubermap.mutants.only %>% 
  ggplot(., mapping = aes(x = score)) + 
  geom_density() + ggtitle("Distribution of E-MAP scores in Gsp1 mutants")

#### look at the distribution of scores for Gsp1 mutants
ubermap.mutants.only %>% 
  ggplot(ubermap.mutants.only, 
         mapping = aes(x = score, 
             y = reorder(Gene_uniq, score, FUN = function(x) quantile(x, prob = 0.05, na.rm = T)))) +
         geom_boxplot() +
        ylab("Mutants ordered by the 5th E-MAP score percentile") +
        xlab("E-MAP score")

ubermap.per.mutant.95conf <- ubermap.mutants.only %>%
  group_by(Gene_uniq) %>%
  summarize(lower_quant = quantile(score, probs = .025, na.rm = T), 
            upper_quant = quantile(score, probs = .975, na.rm = T),
            upper_whisker = boxplot(score, plot = F)$stats[5, ],
            lower_whisker = boxplot(score, plot = F)$stats[1, ]) %>%
  mutate(category = factor(Gene_uniq, levels = Gene_uniq)) %>%
  arrange(lower_quant)
ubermap.per.mutant.95conf_strong <- ubermap.per.mutant.95conf %>% 
  filter(lower_quant < -3)
ubermap.per.mutant.95conf %>% 
  ggplot(., mapping = aes(x = lower_quant, y = upper_quant)) + geom_point() +
  geom_label_repel(data = ubermap.per.mutant.95conf_strong, aes(label = Gene_uniq)) +
  xlab("E-MAP score - lower 2.5 %") + ylab("E-MAP score - higher 97.5 %")
  
ubermap.per.mutant.95conf %>% 
  ggplot(., mapping = aes(x = lower_whisker, y = upper_whisker)) + geom_point() +
  geom_label_repel(data = ubermap.per.mutant.95conf_strong, aes(label = Gene_uniq)) +
  xlab("E-MAP score - lower whisker") + ylab("E-MAP score - upper whisker")


#ubermap.mutants.only <- ubermap.mutants.only %>% 
 # mutate(category = factor(Gene, levels = ubermap.per.mutant.95conf$Gene)) %>%  #### turn into factors so I can order by the lower quartile
#  arrange(category)

## tidy up all but mutants
ubermap.ubergenes.only <- ubermap %>% filter(! (Gene %in% mutants) ) %>%
  #### this is to make sure things like YLR337C - DELTA YLR338W - YLR337C - DELTA YLR338W are not collapse to just YLR337C in the next step
  mutate("Gene" = gsub(Gene, pattern = " - DELTA", replacement = " -- DELTA", perl = T)) %>%
  # this replacement is to remove the redundant ORF before the unique ORF (e.g. YFL008W_TSQ321 - YFL008W_TSQ321 -> keep only the first one) - > first one has the unique i_ for the few redundant ones!
  mutate("Gene" = gsub(Gene, pattern = " - .+$", replacement = "", perl = T))
  ### this is to replace the YOR298C.A with YOR298C-A (R puts . instead of - when loading library genes as column names) - deprecated with tidyverse
  #mutate("Gene" = gsub(Gene, pattern = "\\.", replacement = "-", perl = T))  



(deltaORFS <- ubermap.ubergenes.only %>% 
  filter(., grepl("DELTA", Gene)) %>% 
  pull(Gene) %>% unique())
### most of the DELTA ORFs are dubious reading frames
### I will just use the first ORF for all of them

ubermap.ubergenes.only %>% 
  filter(Gene == "YDR134C -- DELTA YDR133C") %>% 
  mutate("ORF" = gsub(Gene, pattern = "^[0-9]+_", replacement = "", perl = T))
  
#### in this step I make a column that has only the basic ORF - this is necessary so I can match to SGD gene name
ubermap.ubergenes.only <- ubermap.ubergenes.only %>% 
  rename("Gene_uniq" = "Gene") %>%
  mutate("ORF" = gsub(Gene_uniq, pattern = "^[0-9]+_", replacement = "", perl = T)) %>%
  mutate("ORF" = gsub(ORF, pattern = "_.+$", replacement = "", perl = T)) %>% 
  mutate("ORF" = gsub(ORF, pattern = " -- DELTA .+", replacement = "", perl = T))

ubermap.ubergenes.only %>% 
  filter(ORF == "YHR083W")

#### matching to SGD gene names
ubergenes.ubermap.gene_names <- inner_join(ubermap.ubergenes.only, orf_index, by = c("ORF" = "orf")) %>% 
  rename("Gene.gene_name" = "gene_name")
ubergenes.ubermap.gene_names_both <- inner_join(ubergenes.ubermap.gene_names, orf_index, by = c("library.ORF" = "orf")) %>%
  rename("library.gene_name" = "gene_name")
ubergenes.ubermap.gene_names_both %>% 
  filter(ORF == "YKL119C")
write_tsv(ubergenes.ubermap.gene_names_both, path = "basic_E-MAP_data/20180920_preprocessed_ubermap_500_overlap_ubergenes_only.txt")

########### 95% confidence for WT (GSP1-NAT)
(wt <- ubermap.per.mutant.95conf %>% 
    filter(Gene_uniq == "GSP1-NAT") %>% 
    select(-Gene_uniq, -category))
#### 95 % confidence for all the mutants data
(average.95conf <- ubermap.mutants.only %>% 
    summarize(lower_quant = quantile(score, probs = .025, na.rm = T), 
              upper_quant = quantile(score, probs = .975, na.rm = T),
              upper_whisker = boxplot(score, plot = F)$stats[5, ],
              lower_whisker = boxplot(score, plot = F)$stats[1, ]))
(s.lim.point<-c(average.95conf$lower_quant, average.95conf$upper_quant))
plot <- ggplot(ubermap.mutants.only, aes(x = score)) +
  facet_wrap( ~ category) +
  geom_density() +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = lower_quant), size = 1) +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = upper_quant), size = 1) +
  geom_vline(data = wt, aes(xintercept = lower_quant), color = "red", size = 1) + 
  geom_vline(data = wt, aes(xintercept = upper_quant), color = "red", size = 1) +
  geom_vline(data = average.95conf, aes(xintercept = lower_quant), color = "green", size = 1) +
  geom_vline(data = average.95conf, aes(xintercept = upper_quant), color = "green", size = 1) +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = lower_whisker), color = "yellow", size = 1) +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = upper_whisker), color = "yellow", size = 1) +
  geom_vline(data = wt, aes(xintercept = lower_whisker), color = "orange", size = 1) + 
  geom_vline(data = wt, aes(xintercept = upper_whisker), color = "orange", size = 1) +
  geom_vline(data = average.95conf, aes(xintercept = lower_whisker), color = "blue", size = 1) +
  geom_vline(data = average.95conf, aes(xintercept = upper_whisker), color = "blue", size = 1)
print(plot)
pdf(file = "reverse_corr_of_corr/Gsp1_mutants_score_distributions.pdf", width = 14)
print(plot)
dev.off()

####### What is the distribution of scores for the WT construct?
ubermap.mutants.only %>% filter(Gene_uniq == "GSP1-NAT") %>%
  ggplot(aes(x = score)) + geom_histogram() + ggtitle("Distribution of scores for GSP1-NAT")

all.library.genes <- ubermap.mutants.only %>%
  select(library.ORF) %>% unique() %>% pull()
(length(all.library.genes))
##### What is the distribution of scores for the library genes that have a score of lower than 10 in at least one mutant?
highest_mut_emap_score_lib_genes <- ubermap.mutants.only %>% 
  filter(score < -10) %>% 
  pull(library.gene_name) %>% unique()
length(highest_mut_emap_score_lib_genes)
ubermap.mutants.only %>% 
  filter(library.gene_name %in% highest_mut_emap_score_lib_genes) %>%
  ggplot(aes(x = score)) +
  facet_wrap(~library.gene_name) +
  geom_histogram() + ggtitle("Distribution of scores for most negative library genes")

#### get library genes with the most broad score distributions (for Gsp1 mutants genetic interactions)
mut_lib_genes_score_range <- ubermap.mutants.only %>% 
  group_by(library.ORF) %>% 
  do(setNames(data.frame(t(range(.$score, na.rm = T))), c("score_min", "score_max"))) %>% 
  mutate("score_range" = score_max - score_min) %>% 
  arrange(desc(score_range))
mut_lib_genes_score_range %>% 
  ggplot(., aes(x = score_range)) + geom_density()

### there are 304 library gene with a score range larger of 7.5
(mut_lib_genes_score_range_7.5 <- mut_lib_genes_score_range %>% 
  filter(score_range > 7.5) %>% 
  inner_join(., orf_index, by = c("library.ORF" = "orf")) %>% 
  rename("library.gene_name" = "gene_name"))

write_tsv(mut_lib_genes_score_range_7.5, path = "library_genes_with_score_range_over_7.5_in_Gsp1_mut_screens.txt")
######### SIGNIFICANT LIBRARY GENES are defined as as genes that have an E-MAP score outside of the s.lim.point 
## with at least one of the Gsp1 point mutants:
# #### make a subset of ubergenes.ubermap.gene_names_both that only has 
## library genes that have scores outside the s.lim.point with at least one of the mutants
significant.library.genes <- ubermap.mutants.only %>% 
  filter(findInterval(score, s.lim.point) != 1) %>% 
  select(library.ORF) %>% unique() %>% pull()
(length(significant.library.genes))
########### QUESTION: Should I remove from the list of significant.library.genes those that have significant scores with GSP1-NAT (wt control) 
### One could assume that those genetic interaction are an artifact of the construct
(sig.interactions.in.wt <- ubermap.mutants.only %>%
  filter(Gene_uniq == "GSP1-NAT" & findInterval(score, s.lim.point) != 1) %>% 
    inner_join(., orf_index, by = c("library.ORF" = "orf")) %>% 
    inner_join(., mut_lib_genes_score_range, by = "library.ORF"))
###### there is only 7 library genes that are outside of s.lim.point for the WT construct, and 3 of those are nuclear pore proteins
#### they all also have very large score distributions in mutants
###### I would assume it's not a good idea to discard those
##### how does the distribution of scores look like for those library genes ?
ubermap.mutants.only %>%
  filter(library.gene_name %in% sig.interactions.in.wt$library.gene_name) %>%
  ggplot(aes(x = score)) +
  facet_wrap(~library.gene_name) + geom_histogram() + ggtitle("Distribution of E-MAP scores")
###### NOTE: SEM1, SAC3, MFT1 and THP2 are part of THO/TREX complexes
##### Is this relevant?
#### In any case, those 7 library genes all have very broad distribution of scores in mutants, 
#### so don't discard any of them

(ubergenes.ubermap.gene_names_both_sig <- ubergenes.ubermap.gene_names_both %>%
  filter(library.ORF %in% significant.library.genes))
write_tsv(ubergenes.ubermap.gene_names_both_sig, path = "basic_E-MAP_data/preprocessed_ubermap_500_overlap_ubergenes_only_significant.txt")

#### merge the ubergenes.ubermap.gene_names_both and ubermap.mutants.only into one tibble

(combined.emap.data <- bind_rows(ubermap.mutants.only, ubergenes.ubermap.gene_names_both))
tail(combined.emap.data)
write_tsv(combined.emap.data, path = "basic_E-MAP_data/20180921_preprocessed_ubermap_500_overlap_all.txt")
###### only significant subset
(combined.emap.data_sig <- combined.emap.data %>%
    filter(library.ORF %in% significant.library.genes))
write_tsv(combined.emap.data_sig, path = "basic_E-MAP_data/20180921_preprocessed_ubermap_500_overlap_all_sig.txt")

######### non merged ubermap
not.merged.ubermap <- read_tsv("basic_E-MAP_data/ORF_names_ubermap_not_merged.txt",  col_names = T)
#### warning message when reading in the non-merged ubermap
# Warning message:
#   Duplicated column names deduplicated: 'YDR113C' => 'YDR113C_1' [911], 'YLR350W' => 'YLR350W_1' [991], 'YDR432W' => 'YDR432W_1' [1216], 'YGR108W' => 'YGR108W_1' [1303], 'YDR113C' => 'YDR113C_2' [1312], 'YLR350W' => 'YLR350W_2' [1347], 'YGR109C' => 'YGR109C_1' [2039], 'YDR432W' => 'YDR432W_2' [2523], 'YDR113C' => 'YDR113C_3' [2725], 'YER008C' => 'YER008C_1' [2736], 'YPR120C' => 'YPR120C_1' [2823], 'YGR038W' => 'YGR038W_1' [3084], 'YBR160W' => 'YBR160W_1' [3308], 'YBR088C' => 'YBR088C_1' [3371], 'YPR119W' => 'YPR119W_1' [3503], 'YDL155W' => 'YDL155W_1' [3541], 'YDL084W' => 'YDL084W_1' [3548], 'YDR432W' => 'YDR432W_3' [3567], 'YGR038W' => 'YGR038W_2' [3737], 'YJL085W' => 'YJL085W_1' [4140], 'YBR088C' => 'YBR088C_2' [4448], 'YPL153C' => 'YPL153C_1' [4596] 
### tidyverse sorts out duplicated library genes successfully 
not.merged.ubermap %>%
  ggplot(aes(x = YDR113C_1, y = YDR113C)) + geom_point()

#### can't gather before sorting out the Genes that occur multiple times (add a unique identifier to them)
### first find those Genes that are not unique (probably data from separate experiments/labs for the same thing? Stll better not to mix them up)
(not_uniq <- not.merged.ubermap %>% group_by(Gene) %>% 
    summarize("count" = n()) %>% filter(count > 1))

for (temp.gene in not_uniq$Gene) {
  row.names <- which(not.merged.ubermap$Gene == temp.gene)
  for (i in seq_along(row.names)) {
    not.merged.ubermap$Gene[row.names[i]] <- paste(i, not.merged.ubermap$Gene[row.names[i]], sep = "_")
  }
}
### now those few non-unique Genes look like this: 1_YDL155W - YDL155W
#### check
filter(not.merged.ubermap, Gene == "1_YJL085W" | Gene == "2_YJL085W")

# not_uniq is now empty
(not_uniq <- not.merged.ubermap %>% group_by(Gene) %>% 
    summarize("count" = n()) %>% filter(count > 1))
not.merged.ubermap <- not.merged.ubermap %>% 
  gather(library.ORF, score, -Gene) %>%
  drop_na(score) %>%
  rename("query_uniq" = "Gene") %>%
  mutate("query.ORF" = gsub(query_uniq, pattern = " - DELTA", replacement = " -- DELTA", perl = T)) %>% ##### this is to make sure things like YLR337C - DELTA YLR338W - YLR337C - DELTA YLR338W are not collapse to just YLR337C in the next step
  mutate("query.ORF" = gsub(query.ORF, pattern = "^[0-9]+_", replacement = "", perl = T)) %>%
  mutate("query.ORF" = gsub(query.ORF, pattern = "_.+$", replacement = "", perl = T))
filter(not.merged.ubermap, query_uniq != query.ORF)
not.merged.ubermap.gene.names <- inner_join(not.merged.ubermap, orf_index, by = c("query.ORF" = "orf")) %>%
  rename("query_gene_name" = "gene_name")

#### confirm that all query_uniq are unique
##########################
# ubermap_per_lib_orf_count <- not.merged.ubermap.gene.names %>% group_by(library.ORF) %>%
#   summarise("gene_count" = n())
# unique(ubermap_per_lib_orf_count$gene_count)
# ubermap_per_gene_count <- not.merged.ubermap.gene.names %>% group_by(query_uniq) %>%
#   summarise("lib_count" = n())
# unique(ubermap_per_gene_count$lib_count)
# ubermap_per_gene_count %>% 
#   filter(lib_count > mean(ubermap_per_gene_count$lib_count)) %>%
#   pull(query_uniq)
# ubermap_per_gene_count %>% filter(lib_count == 3072)
save(not.merged.ubermap.gene.names, file = "basic_E-MAP_data/ubermap_not_merged.RData")
##########################################

#### preprocess histone point mutations data
histones <- read_delim("basic_E-MAP_data/20180226_histones_pE-MAPs.txt", delim = "\t", col_names = T)
histones <- histones %>% gather(library.gene_name, score, -Gene)
write_tsv(histones, path = "basic_E-MAP_data/histones_pE-MAPs.txt")
