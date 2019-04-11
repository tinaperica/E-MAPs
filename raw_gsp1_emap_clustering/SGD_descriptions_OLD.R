#
# title: "get_SGD_descriptions"
# author: "Chris Mathy"
# date: "4/2/2019"
# output: html_document

#This script creates a data table containing the descriptions
#from the SGD for each gene, labeled by its name and ORF.

# BiocManager::install("org.Sc.sgd.db", version = "3.8")
library(org.Sc.sgd.db)


# Get descriptions

description_map <- org.Sc.sgdDESCRIPTION
ORFs <- mappedkeys(description_map) 
descriptions <- as.list(description_map[ORFs])
description_df <- tibble(ORF = names(descriptions), description = unlist(descriptions))


# Get common names

names_map <- org.Sc.sgdGENENAME
ORFs <- mappedkeys(names_map) 
names <- as.list(names_map[ORFs])
names_df <- tibble(ORF = names(names), name = unlist(names))


# Join and save

SGD_descriptions <-
  inner_join(description_df, names_df, by = 'ORF') %>% 
  dplyr::select(name, ORF, description)

write_delim(SGD_descriptions, './SGD_descriptions_from_SGD_Bioconductor.txt', delim='\t')


