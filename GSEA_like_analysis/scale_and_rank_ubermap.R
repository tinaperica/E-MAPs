library(tidyverse)
load("basic_E-MAP_data/20190129_ubermap_for_correlations.RData")

ubermap %>% ggplot(aes(x = score)) + geom_density()
Gsp1_emap <- ubermap %>% 
  filter(query_ORF == "YLR293C" & query_uniq != "GSP1_tsa256")
sga <- ubermap %>% 
  filter(query_ORF != "YLR293C" | query_uniq == "GSP1_tsa256")
scaled_pos_Gsp1_emap <- Gsp1_emap %>% 
  filter(score > 0) %>% 
  mutate("scaled_score" = scale(score))
scaled_neg_Gsp1_emap <- Gsp1_emap %>% 
  filter(score < 0) %>% 
  mutate("scaled_score" = scale(abs(score)))
scaled_pos_sga <- sga %>% 
  filter(score > 0) %>% 
  mutate("scaled_score" = scale(score))
scaled_neg_sga <- sga %>% 
  filter(score < 0) %>% 
  mutate("scaled_score" = scale(abs(score)))
scaled_ubermap <- bind_rows(scaled_pos_Gsp1_emap, scaled_neg_Gsp1_emap, scaled_pos_sga, scaled_neg_sga) %>% 
  arrange(desc(scaled_score))
scaled_ubermap %>% ggplot(aes(x = score)) + geom_density()

save(scaled_ubermap, file = "basic_E-MAP_data/20190129_scaled_sga_ubermap_for_GSEA.RData")
