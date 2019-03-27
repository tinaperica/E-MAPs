library(tidyverse)

load("basic_E-MAP_data/spitzemapko_for_corr.rda")


spitzemapko_for_gia <- spitzemapko_for_corr %>% 
  rename("Gene" = query_allele_name) %>% 
  select(-weight) %>% 
  spread(array_ORF, score)

spitzemapko_for_gia[spitzemapko_for_gia > -4 & spitzemapko_for_gia < 4] <- 0

columns_to_drop <- spitzemapko_for_gia %>% 
  gather(library_gene, score, -Gene) %>% 
  group_by(library_gene) %>% 
  summarise("sum_score" = sum(abs(score), na.rm = T)) %>% 
  arrange(sum_score) %>% 
  filter(sum_score == 0 | is.na(sum_score)) %>% 
  pull(library_gene) %>% unique()
#### there are no columns to drop

write_tsv(spitzemapko_for_gia, "gia_analysis/spitzemapko_for_gia.txt", na = "")

#### left annotation for gia
gsp1_left_anno <- spitzemapko_for_gia %>% 
  filter(grepl("GSP1", Gene)) %>% 
  select(Gene) %>% 
  unique() %>% 
  separate(col = Gene, into = c("gene", "mutant"), sep = " - ", remove = F) %>% 
  select("X1" = Gene, "X2" = mutant)
sga_left_anno <- spitzemapko_for_gia %>% 
  filter(! grepl("GSP1", Gene)) %>% 
  select(Gene) %>% 
  unique() %>% 
  mutate("X2" = str_to_upper(Gene)) %>% 
  rename("X1" = Gene)
gia_left_annotation <- bind_rows(gsp1_left_anno, sga_left_anno)
write_tsv(gia_left_annotation, "gia_analysis/spitzemapko_left_annotation_for_gia.txt")
