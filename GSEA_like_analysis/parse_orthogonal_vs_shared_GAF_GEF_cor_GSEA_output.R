#### for parsing and analysing the output from GSEA 3.0
library(tidyverse)
#### orthogonal positive is quadrant IV
### orthogonal negative is quadrant II
### shared positive is quadrant I
## shared negative is quadrant III
input <- "20190321"
input_dir <- file.path("GSEA_like_analysis/SGSEA 3.0 output folder", input)

gsea_out_folders <- dir(input_dir)
files_tib <- tibble("diagonal" = as.character(), "filepaths" = as.character(), "enrichment_sign" = as.character())
for (o in seq_along(gsea_out_folders)) {
  diagonal <- gsea_out_folders[o] %>% str_split("_") %>% .[[1]]  %>% .[3]
  dir_path <- file.path(input_dir, gsea_out_folders[o])
  split_dir <- str_split(gsea_out_folders[o], pattern = ".GseaPreranked.")
  filenames <- dir(dir_path)
  filenames <- filenames[grepl(filenames, pattern = "gsea_report_for_na.+xls", perl = T)]
  enrichment_sign <- filenames %>% str_split("_") %>% unlist() %>% .[c(5, 11)]
  filepaths <- file.path(dir_path, filenames)
  files_tib <- bind_rows(files_tib, tibble("diagonal" = diagonal, "filepaths" = filepaths, 
                                           "enrichment_sign" = enrichment_sign))
}
enrichment_data <- tibble("gene_set" = as.character(), "gene_set_size"= integer(), "norm_ES" = double(), 
                          "FDR" = double(), "diagonal" = character(), "sign" = character())

for (i in seq_along(files_tib$filepaths))  {
  filepath <- files_tib$filepaths[i]
  diagonal <- files_tib$diagonal[i]
  enrichment_sign <- files_tib$enrichment_sign[i]
  input_tib <- read_tsv(filepath, col_names = T) %>% 
    select("gene_set" = NAME, "gene_set_size" = SIZE, 'norm_ES' = NES, "FDR" = `FDR q-val`) %>% 
    mutate("diagonal" = diagonal, "sign" = enrichment_sign)
  enrichment_data <- bind_rows(enrichment_data, input_tib) %>% 
    arrange(gene_set)
}

enrichment_data <- enrichment_data %>% 
  mutate("quadrant" = ifelse( (diagonal == 'orthogonal' & sign == 'pos'), 'IV', NA)) %>% 
  mutate("quadrant" = ifelse( (diagonal == 'orthogonal' & sign == 'neg'), 'II', quadrant)) %>% 
  mutate("quadrant" = ifelse( (diagonal == "shared" & sign == "pos"), "I", quadrant )) %>% 
  mutate("quadrant" = ifelse( (diagonal == "shared" & sign == "neg"), "III", quadrant ))

to_save <- enrichment_data %>% filter(FDR < 0.1) %>% 
  arrange(quadrant, FDR)
write_tsv(to_save, "GSEA_like_analysis/enrichments_by_quadrant_of_GAP_GEF_rank_correlation_of_library_genes.txt")
