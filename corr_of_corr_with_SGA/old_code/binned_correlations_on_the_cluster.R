library(tidyverse)
load("basic_E-MAP_data/20181230_ubermap_for_correlations.RData")
head(ubermap)
query_uniq_to_select <- c("T34A", "R108L", "T34E", "T34G", "D79S", "H141R", "GSP1-NAT", "T137G", "R78K",
                          "Y157A", "K101R", "YLR293C_tsa256", "Q147E", "R112S", "G80A")
core_partner_gene_names <- c("SRM1", "MOG1", "YRB1","YRB2", "YRB30", "RNA1", "LOS1", "MSN5", "MTR10",
                             "CSE1", "CRM1", "CEX1", "SRP1", "KAP95", "KAP120", "KAP123", "KAP104", 
                             "KAP114", "PSE1", "NMD5", "NUP42", "NTF2")
test_ubermap <- ubermap %>% 
  filter(query_uniq %in% query_uniq_to_select | query_gene_name %in% core_partner_gene_names)

