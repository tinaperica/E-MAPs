library(tidyverse)

load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
e.map <- as_tibble(e.map)
e.map %>% pull(mutant) %>% unique()
