# extra code
# 


# this code shows that there are no duplicate measurements in the gsp1 data (after averaging, as expected)
# e.map %>% select(mutant, library_gene_name) %>% unique() %>% select(mutant) %>%  unlist() %>% length()

gsp1_spread <-
  e.map %>%
  as.tibble %>% 
  select(-library_ORF) %>% 
  spread(library_gene_name, score)


ptm <- proc.time()
proc.time() - ptm
