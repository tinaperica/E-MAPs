library(tidyverse)
### SGA nomenclature
## _tsa - TS mutant in array (library, Kan) strain
## _dma - gene deletion in array (library, Kan) strain
## _tsg - TS mutant in query (query, NAT) strain
## _sn - gene deletion in query (query, NAT) strain
load("cF3_scaled_SGA2016.RData")
load("June2016_Gsp1_E-MAP_data.RData")
ORF_gene_SGD <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F) %>%
  select("ORF" = X1, "gene_name"= X2) %>% unique()
  
e.map <- as_tibble(e.map) %>% 
  rename("query_uniq" = mutant) %>% 
  mutate("library_uniq" = library_ORF) %>% 
  mutate("query_ORF" = "YLR293C", "query_gene_name" = "GSP1")
  
e.map.library_ORFs <- e.map %>% pull(library_ORF) %>% unique()
SGA_scaled <- SGA_scaled %>% 
  select(Query_Strain_ID, Query_ORF, Array_Strain_ID, Array_ORF, sga_score_scaled) %>% 
  filter(Query_ORF %in% e.map.library_ORFs | Array_ORF %in% e.map.library_ORFs)

flipped_SGA_scaled <- SGA_scaled %>% 
  select( Query_Strain_ID = Array_Strain_ID, Query_ORF = Array_ORF,
          Array_Strain_ID = Query_Strain_ID, Array_ORF = Query_ORF,
          sga_score_scaled = sga_score_scaled)
SGA_scaled <- SGA_scaled %>% 
  bind_rows(., flipped_SGA_scaled) %>% 
  filter(Array_ORF %in% e.map.library_ORFs) %>% 
  filter( ! grepl("ts", Array_Strain_ID)) %>% unique() %>% 
  inner_join(., ORF_gene_SGD, by = c("Query_ORF" = "ORF")) %>% 
  rename("query_gene_name" = gene_name,
         "query_uniq" = Query_Strain_ID, "query_ORF" = Query_ORF) %>% 
  inner_join(., ORF_gene_SGD, by = c("Array_ORF" = "ORF")) %>% 
  rename("library_gene_name" = gene_name, "library_ORF" = Array_ORF, 
         "library_uniq" = Array_Strain_ID,
         "score" = sga_score_scaled)

ubermap <- full_join(SGA_scaled, e.map)
save(ubermap, file = "SGA_Gsp1emap_merged.RData")
