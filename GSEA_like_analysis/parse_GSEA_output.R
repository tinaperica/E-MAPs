#### for parsing and analysing the output from GSEA 3.0
library(tidyverse)

input <- "20190307"
input_dir <- file.path("GSEA_like_analysis/SGSEA 3.0 output folder", input)

gsea_out_folders <- dir(input_dir)

for (o in seq_along(gsea_out_folders)) {
  dir_path <- file.path(input_dir, gsea_out_folders[o])
  split_dir <- str_split(gsea_out_folders[o], pattern = ".GseaPreranked.")
  files <- dir(dir_path)
  files <- files[grepl(files, pattern = "gsea_report_for_na.+xls", perl = T)]
}
for (i in seq_along(files))  {    ## this loops through files (one per cluster)
  filepath <- file.path(corr_of_corr_path, files[i])
  load(filepath)
  corr_of_corr <- as_tibble(corr_of_corr) %>% 
    filter(corr_of_corr_type == "w_pearson_corr") %>% 
    mutate("fdr" = p.adjust(p.value, method = "fdr")) %>% 
    #filter(fdr < 0.05) %>% 
    inner_join(., Gene_uniq_matcher, by = c("geneB" = "Gene_uniq")) %>% 
    select(geneA, geneB, Gene.gene_name, cluster, corr, fdr)
  combined_data <- bind_rows(combined_data, corr_of_corr)
}