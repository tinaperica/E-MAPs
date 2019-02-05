library(tidyverse)
load("basic_E-MAP_data/20181230_ubermap_for_correlations.RData")
head(ubermap)
query_uniq_to_select <- c("T34A", "R108L", "T34E", "T34G", "D79S", "H141R", "GSP1-NAT", "T137G", "R78K",
                          "Y157A", "K101R", "YLR293C_tsa256", "Q147E")
test_core_partner_gene_names <- c("SRM1", "MOG1", "YRB1", "RNA1", "MSN5", "MTR10",
                             "CSE1", "CRM1", "CEX1", "SRP1", "KAP95", 
                             "KAP114", "PSE1", "NUP42", "NTF2")
test_ubermap <- ubermap %>% 
  filter(query_uniq %in% query_uniq_to_select | query_gene_name %in% test_core_partner_gene_names) %>% 
  select(query_uniq, library_ORF, score)

test_query_uniq <- test_ubermap %>% pull(query_uniq) %>% unique()
all_test_pairs <- combn(test_query_uniq, 2)
n_pairs <- length(all_test_pairs)/2

for (i in 1:n_pairs ) {
  query1 <- all_test_pairs[1, i]
  query1_emap <- test_ubermap %>% filter(query_uniq == query1)
  query2 <- all_test_pairs[2, i]
  query2_emap <- test_ubermap %>% filter(query_uniq == query2)
  pair_emap <- query1_emap %>% 
    inner_join(., query2_emap, by = "library_ORF" ) %>% 
    filter(score.x > 1 | score.y > 1 | score.x < -2 | score.y < -2 )
  pair_emap %>% ggplot(aes(x = score.x, y = score.y)) + 
    geom_point() + xlab(query1) + ylab(query2) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
  
  quadrants <- as.tibble(table(pair_emap$score.x > 0, pair_emap$score.y > 0)) %>% 
    mutate("quadrant" = NA) %>% 
    mutate("quadrant" = ifelse((Var1 == "FALSE" & Var2 == "FALSE"), 3, quadrant)) %>% 
    mutate("quadrant" = ifelse((Var1 == "TRUE" & Var2 == "TRUE"), 1, quadrant)) %>%
    mutate("quadrant" = ifelse((Var1 == "FALSE" & Var2 == "TRUE"), 2, quadrant)) %>%
    mutate("quadrant" = ifelse((Var1 == "TRUE" & Var2 == "FALSE"), 4, quadrant))
    
  
}