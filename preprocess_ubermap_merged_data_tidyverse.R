library(tidyverse)
library(gmp)
bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')
##### bad strains are library strains to be removed (according to Amelie's file, and Hannes's instruction, those should be removed)

### merged and averaged E-MAP data for Gsp1 mutants
#### read it in and make it tidy
e.map <- read_delim("basic_E-MAP_data/avg_merged_June2016_screen_for_Gia.txt", col_names = T, delim = "\t")   ## use export for Gia because it has ORF names for library
e.map <- e.map %>% gather(library, score, -Gene) %>%
  filter(! library %in% bad_strains) %>%
  arrange(desc(Gene), desc(library)) %>%
  mutate("Gene" = gsub(Gene, pattern = "GSP1:", replacement = "", perl = T))

Gsp1_queries <- unique(e.map[["Gene"]])
### ORF to gene name annotation from SGD
orf_gene_name_index <- read_delim("orf_gene_GO_sgd_annotation.txt", col_names = F, delim = "\t")
orf_index <- unique(tibble("orf" = orf_gene_name_index$X1, "gene_name" = orf_gene_name_index$X2))
rm(orf_gene_name_index)
#############################

### Gsp1 mutants merged with ubermap (gene names) - this is only necessary to distinguish point mutants (they are all the same in the ORF file)
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
### first find those Genes that are not unique (probably data from separate experiments/labs for the same thing? Stll better not to mix them up)
(not_uniq <- ubermap %>% group_by(Gene) %>% 
  summarize("count" = n()) %>% filter(count > 1))

for (temp.gene in not_uniq$Gene) {
  row.names <- which(ubermap$Gene == temp.gene)
  for (i in seq_along(row.names)) {
    ubermap$Gene[row.names[i]] <- paste(i, ubermap$Gene[row.names[i]], sep = "_")
  }
}
### now those few non-unique Genes look like this: 1_YDL155W - YDL155W
#### check
filter(ubermap, Gene == "1_YDL155W - YDL155W" | Gene == "2_YDL155W - YDL155W")
(not_uniq <- ubermap %>% group_by(Gene) %>% 
    summarize("count" = n()) %>% filter(count > 1))

ubermap <- ubermap %>% gather(library.ORF, score, -Gene)

## mutants only
ubermap.mutants.only <- ubermap %>% filter( Gene %in% mutants) %>%
  mutate("library.ORF" = gsub(library.ORF, pattern = "\\.", replacement = "-", perl = T)) %>% ### this is to replace the YOR298C.A with YOR298C-A (R puts . instead of - when loading library genes as column names)
  add_column("ORF" = "YLR293C")
write_delim(ubermap.mutants.only, path = "basic_E-MAP_data/preprocessed_ubermap_mut_only.txt", delim = "\t")
ggplot(ubermap.mutants.only, mapping = aes(x = score)) + geom_density()
ubermap.mutants.only %>% 
  ggplot(ubermap.mutants.only, mapping = aes(x = score, y = reorder(Gene, score, FUN = function(x) quantile(x, prob = 0.025, na.rm = T)))) +
         geom_boxplot()
ubermap.per.mutant.95conf <- ubermap.mutants.only %>%
  group_by(Gene) %>%
  summarize(lower = quantile(score, probs = .025, na.rm = T), 
            upper = quantile(score, probs = .975, na.rm = T)) %>%
  mutate(category = factor(Gene, levels = Gene)) %>%
  arrange(lower)
ubermap.mutants.only <- ubermap.mutants.only %>% 
  mutate(category = factor(Gene, levels = ubermap.per.mutant.95conf$Gene)) %>%  #### turn into factors so I can order by the lower quartile
  arrange(category)

## tidy up all but mutants
ubermap.ubergenes.only <- ubermap %>% filter(! (Gene %in% mutants) ) %>%
  mutate("Gene" = gsub(Gene, pattern = " - DELTA", replacement = " -- DELTA", perl = T)) %>% ##### this is to make sure things like YLR337C - DELTA YLR338W - YLR337C - DELTA YLR338W are not collapse to just YLR337C in the next step
  mutate("Gene" = gsub(Gene, pattern = " - .+$", replacement = "", perl = T)) %>% ### this replacement is to remove the redundant ORF before the unique ORF (e.g. YFL008W_TSQ321 - YFL008W_TSQ321 -> keep only the first one) - > first one has the unique i_ for the few redundant ones!
  mutate("Gene" = gsub(Gene, pattern = "\\.", replacement = "-", perl = T))  ### this is to replace the YOR298C.A with YOR298C-A (R puts . instead of - when loading library genes as column names)
rm(ubermap)
filter(ubermap.ubergenes.only, grepl("DELTA", Gene))
#### in this step I make a column that has only the basic ORF - this is necessary so I can match to SGD gene name
ubermap.ubergenes.only <- ubermap.ubergenes.only %>% rename("Gene_uniq" = "Gene") %>%
    mutate("ORF" = gsub(Gene_uniq, pattern = "^[0-9]+_", replacement = "", perl = T)) %>%
    mutate("ORF" = gsub(ORF, pattern = "_.+$", replacement = "", perl = T))
#### matching to SGD gene names
ubergenes.ubermap.gene_names <- inner_join(ubermap.ubergenes.only, orf_index, by = c("ORF" = "orf")) %>% 
  rename("Gene.gene_name" = "gene_name")
ubergenes.ubermap.gene_names_both <- inner_join(ubergenes.ubermap.gene_names, orf_index, by = c("library.ORF" = "orf")) %>%
  rename("library.gene_name" = "gene_name")
write_delim(ubergenes.ubermap.gene_names_both, path = "basic_E-MAP_data/preprocessed_ubermap_ubergenes_only.txt", delim = "\t")

########### 95% confidence for WT (GSP1-NAT)
(wt <- ubermap.per.mutant.95conf %>% filter(Gene == "GSP1-NAT") %>% select(-Gene, -category))
#### 95 % confidence for all the mutants data
(average.95conf <- ubermap.mutants.only %>% 
    summarize(lower = quantile(score, probs = .025, na.rm = T),
              upper = quantile(score, probs = .975, na.rm = T)))
(s.lim.point<-c(average.95conf$lower, average.95conf$upper))
plot <- ggplot(ubermap.mutants.only, aes(x = score)) +
  facet_wrap( ~ category) +
  geom_density() +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = lower)) +
  geom_vline(data = ubermap.per.mutant.95conf, aes(xintercept = upper)) +
  geom_vline(data = wt, aes(xintercept = lower), color = "red") + 
  geom_vline(data = wt, aes(xintercept = upper), color = "red") +
  geom_vline(data = average.95conf, aes(xintercept = lower), color = "green") +
  geom_vline(data = average.95conf, aes(xintercept = upper), color = "green")
print(plot)
pdf(file = "Gsp1_mutants_score_distributions.pdf", width = 14)
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
length(highest_mut_emap_score_ib_genes)
filter(ubermap.mutants.only, library.gene_name %in% highest_mut_emap_score_lib_genes) %>%
  ggplot(aes(x = score)) +
  facet_wrap(~library.gene_name) +
  geom_histogram() + ggtitle("Distribution of scores for most negative library genes")
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
  filter(Gene_uniq == "GSP1-NAT" & findInterval(score, s.lim.point) != 1))
###### there is only 7 library genes that are outside of s.lim.point for the WT construct, and 3 of those are nuclear pore proteins
###### I would assume it's not a good idea to discard those
##### how do the scores with those library genes look with the mutants?
ubermap.mutants.only %>%
  filter(library.gene_name %in% sig.interactions.in.wt$library.gene_name) %>%
  ggplot(aes(x = score)) +
  facet_wrap(~library.gene_name) + geom_histogram() + ggtitle("Distribution of E-MAP scores")
###### NOTE: SEM1, SAC3, MFT1 and THP2 are part of THO/TREX complexes
##### Is this relevant?
#### In any case, those 7 library genes all have very broad distribution of scores in mutants, so don't discard any of them

(ubergenes.ubermap.gene_names_both_sig <- ubergenes.ubermap.gene_names_both %>%
  filter(library.ORF %in% significant.library.genes))
write_delim(ubergenes.ubermap.gene_names_both_sig, path = "basic_E-MAP_data/preprocessed_ubermap_ubergenes_only_significant.txt", delim = "\t")

#### merge the ubergenes.ubermap.gene_names_both and ubermap.mutants.only into one tibble
ubermap.mutants.only <- ubermap.mutants.only %>%
  select("Gene_uniq" = Gene, library.ORF, score, ORF) %>%
  mutate("Gene.gene_name" = "GSP1")
ubermap.mutants.only <- ubermap.mutants.only %>%
  inner_join(orf_index, by = c("library.ORF" = "orf")) %>%
  rename("library.gene_name" = gene_name)

(combined.emap.data <- bind_rows(ubermap.mutants.only, ubergenes.ubermap.gene_names_both))
tail(combined.emap.data)
write_delim(combined.emap.data, path = "basic_E-MAP_data/preprocessed_ubermap_all.txt", delim = "\t")
###### only significant subset
(combined.emap.data_sig <- combined.emap.data %>%
    filter(library.ORF %in% significant.library.genes))
write_delim(combined.emap.data_sig, path = "basic_E-MAP_data/preprocessed_ubermap_all_sig.txt", delim = "\t")

######### non merged ubermap
not.merged.ubermap <- read_delim("basic_E-MAP_data/ORF_names_ubermap_not_merged.txt", 
                                 delim = "\t", col_names = T)
#### warning message when reading in the non-merged ubermap
# Warning message:
#   Duplicated column names deduplicated: 'YDR113C' => 'YDR113C_1' [911], 'YLR350W' => 'YLR350W_1' [991], 'YDR432W' => 'YDR432W_1' [1216], 'YGR108W' => 'YGR108W_1' [1303], 'YDR113C' => 'YDR113C_2' [1312], 'YLR350W' => 'YLR350W_2' [1347], 'YGR109C' => 'YGR109C_1' [2039], 'YDR432W' => 'YDR432W_2' [2523], 'YDR113C' => 'YDR113C_3' [2725], 'YER008C' => 'YER008C_1' [2736], 'YPR120C' => 'YPR120C_1' [2823], 'YGR038W' => 'YGR038W_1' [3084], 'YBR160W' => 'YBR160W_1' [3308], 'YBR088C' => 'YBR088C_1' [3371], 'YPR119W' => 'YPR119W_1' [3503], 'YDL155W' => 'YDL155W_1' [3541], 'YDL084W' => 'YDL084W_1' [3548], 'YDR432W' => 'YDR432W_3' [3567], 'YGR038W' => 'YGR038W_2' [3737], 'YJL085W' => 'YJL085W_1' [4140], 'YBR088C' => 'YBR088C_2' [4448], 'YPL153C' => 'YPL153C_1' [4596] 

#### can't gather before sorting out those few Genes that occur multiple times (add a unique identifier to them)
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
write_delim(histones, path = "basic_E-MAP_data/histones_pE-MAPs.txt", delim = "\t")
