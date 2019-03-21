library(tidyverse)

### this code makes gene groups based on Wodak complexes with more than 3 members and the GO slims terms
#### to use for the gene enrichment analyses

### first make clusters based on GO slims
GO_slims <- read_delim( "go_slim_mapping.tab.txt", col_names = F, delim = "\t") %>% 
  select("ORF" = X1, "Gene" = X2, "GO_Slim_term" = X5) %>% 
  mutate("Gene" = ifelse(is.na(Gene), ORF, Gene))

orf_gene_name_index <- GO_slims %>% select(ORF, Gene) %>% unique()
write_tsv(orf_gene_name_index, path = "orf_gene_GO_sgd_annotation.txt", col_names = F)
#### check for overly general slims
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
    summarize("count" = n()) %>% arrange(desc(count)))
### remove some terms
##############################################
GO_slim_terms_to_remove <- as.character(expression(
  signaling, biological_process, response, "response to chemical", "protein complex biogenesis", other, "not_yet_annotated",
  "ion binding", "cellular_component", "molecular_function", membrane, cytoplasm, nucleus, "structural molecule activity", "ATPase activity"
))
GO_slims <- GO_slims %>% filter( ! (GO_Slim_term %in% GO_slim_terms_to_remove))
(GO_slims_count <- GO_slims %>% group_by(GO_Slim_term) %>%
    summarize("count" = n()) %>% arrange(desc(count)))
tail(GO_slims_count)
ggplot(GO_slims_count, aes(x = count)) + geom_histogram() + ggtitle("Go slims terms / cluster sizes")
GO_slim_terms <- GO_slims %>% pull(GO_Slim_term) %>% unique()

## now add Wodak complexes with more than 3 members
complex_tibble <- read_tsv("titration_curves/all_wodak_complexes.txt") %>% 
  rename("Gene" = Name)
complexes <- complex_tibble %>% pull(Complex) %>% unique()
### make each remaining GO_slim term it's own group
#### and then make some composite groups (e.g. combin ER and Golgi, or combine multiple terms into mRNA processing group etc.)
##################################################################

### finally add some manually made gene groups
manual <- read_tsv("GSEA_like_analysis/manual_gene_groups.txt", col_names = F) %>% 
  rename("ORF" = X1, "Gene" = X2, "term" = X3)
manual_terms <- manual %>% pull(term) %>% unique()
## first add all the terms as groups
gene_groups <- list()
for (i in seq_along(GO_slim_terms)) {
  gene_groups[[GO_slim_terms[i]]] <- GO_slim_terms[i]
}
for (i in seq_along(complexes)) {
  gene_groups[[complexes[i]]] <- complexes[i]
}
for (i in seq_along(manual_terms)) {
  gene_groups[[manual_terms[i]]] <- manual_terms[i]
}
### then manually select some terms based on which I will make GO based clusters
#### make sure that composite GO_slims_terms have unique names, that don't already exist as GO_slims_terms
##### this approach gives a small number of relativelly general clusters
gene_groups[["ribosomes and translation"]] <- unique ( GO_slims$GO_Slim_term[ 
  grep( paste( c("rRNA", "ribosom", "translation", "tRNA", "snoRNA"), collapse = "|"), GO_slims$GO_Slim_term)
  ])
gene_groups[["transcription and mRNA processing"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mRNA", "RNA splicing", "transcription", "RNA modification", "mRNA processing",
                 "rRNA processing", "tRNA processing", "snoRNA processing"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["transcription"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("transcription"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["RNA processing"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mRNA processing", "splicing", "snoRNA processing"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["Golgi and ER"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("protein lipidation", "protein maturation", "endocytosis", "regulation of transport", "glycosylation",
                 "vesicle organization","endosom", "Golgi", "endoplasmic"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["Golgi"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("Golgi"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["ER"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("endoplasmic reticulum"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["peroxisomes"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("peroxisom"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["vacuoles"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("vacuol"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["mitochondria"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("mitochond", "respirat"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["chromatin"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("chromatin","histone", "telomere", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["cytoskeleton and microtubules"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cytoskelet", "cytokinesis", "chromosome segregation", "microtubule"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["cell cycle"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("cell cycle", "recombination", "chromosome segregation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["budding"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("bud", "fission", "sporulation"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["lipids"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("lipid"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["nuclear transport and organization"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("nuclear transport", "nucleus"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
gene_groups[["metabolic"]] <- unique ( GO_slims$GO_Slim_term[
  grep( paste( c("metabolic", "metabolite"), collapse = "|"), GO_slims$GO_Slim_term) 
  ])
#####################

(groups <- names(gene_groups))

term_orf_index <- complex_tibble %>% 
  rename("GO_Slim_term" = "Complex") %>% 
  bind_rows(., GO_slims) %>% 
  rename("term" = GO_Slim_term) %>% 
  bind_rows(., manual)
gene_categories <- tibble("ORF" = character(), Gene = "character", term = "character")
for (i in seq_along(groups)) {
  terms <- gene_groups[[groups[i]]]
  temp <- term_orf_index %>% 
    filter(term %in% terms)
  gene_categories <- bind_rows(gene_categories, temp)
}


larger_than_5_terms <- gene_categories %>% 
  group_by(term) %>% 
  summarize("count" = n()) %>% 
  filter(count > 5) %>% 
  pull(term) %>% unique()
gene_categories <- gene_categories %>% 
  filter(term %in% larger_than_5_terms)
write_tsv(gene_categories, path = "GSEA_like_analysis/gene_groups.txt")

for_gia <- gene_categories %>% select(ORF, term)
write_tsv(for_gia, "GSEA_like_analysis/gia_analysis/gia_right_annotations_gene_groups.txt")
### output the categories for GSEA 3.0 Java script (output in .gmt format)
### .gmt format defined in http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#TXT:_Text_file_format_for_expression_dataset_.28.2A.txt.29


gene_categories <- gene_categories %>% 
  mutate("gene_set_type" = ifelse(term %in% complexes, "Wodak complex", "GO slims"))
terms <- gene_categories %>% pull(term) %>% unique()
if (file.exists("GSEA_like_analysis/gene_groups.gmt")) {
  unlink("GSEA_like_analysis/gene_groups.gmt")
}
for (i in seq_along(terms)) {
  gene.set <- gene_categories %>% filter(term == terms[i]) %>% 
    pull(Gene)
  cat <- gene_categories %>% filter(term == terms[i]) %>% pull(gene_set_type) %>% unique()
  gmt.row <- c(terms[i], cat, gene.set)
  write(gmt.row, file = "GSEA_like_analysis/gene_groups.gmt", ncolumns = length(gmt.row), append = T, sep = "\t")
}
