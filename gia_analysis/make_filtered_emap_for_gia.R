#### prepare the raw E-MAP scores for gia by setting all interactions between -2 and 2 to 0
### I only want enrichments based on significant S-scores
library(tidyverse)
bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')

emap <- read_tsv("gia_analysis/avg_merged_June2016_screen_for_Gia.txt")
emap <- emap %>% select(-one_of(bad_strains))

emap[emap > -2 & emap < 2] <- 0
write_tsv(emap, "gia_analysis/emap_screen_filtered.txt", na = "")
