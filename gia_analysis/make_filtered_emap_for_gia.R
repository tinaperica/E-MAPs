#### prepare the raw E-MAP scores for gia by setting all interactions between -2 and 2 to 0
### I only want enrichments based on significant S-scores
library(tidyverse)
bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')

emap <- read_tsv("gia_analysis/avg_merged_June2016_screen_for_Gia.txt")
emap <- emap %>% select(-one_of(bad_strains))

emap[emap > -2 & emap < 2] <- 0
columns_to_drop <- emap %>% 
  gather(library_gene, score, -Gene) %>% 
  group_by(library_gene) %>% 
  summarise("sum_score" = sum(score)) %>% 
  filter(sum_score == 0) %>% 
  pull(library_gene) %>% unique()
emap <- emap %>% select(-one_of(columns_to_drop))
write_tsv(emap, "gia_analysis/emap_screen_filtered_at_50perc_conf.txt", na = "")

emap[emap > -4 & emap < 4] <- 0

columns_to_drop <- emap %>% 
  gather(library_gene, score, -Gene) %>% 
  group_by(library_gene) %>% 
  summarise("sum_score" = sum(score)) %>% 
  filter(sum_score == 0) %>% 
  pull(library_gene) %>% unique()
emap <- emap %>% select(-one_of(columns_to_drop))
write_tsv(emap, "gia_analysis/emap_screen_filtered_at_85perc_conf.txt", na = "")

### GI from SGD
Gsp1_sgd_gi_and_ppi <- read_tsv("~/Desktop/GSP1_interactions.txt") %>% 
  select(`Interactor Systematic Name_1`, Type) %>% 
  pull(`Interactor Systematic Name_1`) %>% unique()
### our GI
our_gi <- emap %>% gather(library_gene, score, -Gene) %>% pull(library_gene) %>% unique()
length(intersect(our_gi, Gsp1_sgd_gi_and_ppi))
