---
title: "costanzo2016_processing"
author: "Chris Mathy"
date: "11/15/2018"
output: html_document
---

# import packages
library(tidyverse)

# download SGA data from costanzo2016

url <- 'http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/Data%20File%20S1_Raw%20genetic%20interaction%20datasets:%20Pair-wise%20interaction%20format.zip'
download_dir <- dir.create('costanzo2016_data/', showWarnings = FALSE)
download.file(url=url, destfile='costanzo2016_data/costanzo2016.zip')
unzip('costanzo2016_data/costanzo2016.zip', exdir='costanzo2016_data')
file.rename(from='costanzo2016_data/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format', to='costanzo2016_data/datafileS1')

# read in SGA data from costanzo 2016

SGA_ExE <- read_delim('costanzo2016_data/datafileS1/SGA_ExE.txt', delim='\t') # essential vs. essential
SGA_ExN <- read_delim('costanzo2016_data/datafileS1/SGA_ExN.txt', delim='\t') # essential vs. non-essential
SGA_NxN <- read_delim('costanzo2016_data/datafileS1/SGA_NxN.txt', delim='\t') # non-essential vs. non-essential
SGA_DAmP <- read_delim('costanzo2016_data/datafileS1/SGA_DAmP.txt', delim='\t') # non-essential vs. non-essential

SGA <-
  bind_rows('ExE' = SGA_ExE,
            'ExN' = SGA_ExN,
            'NxN' = SGA_NxN,
            'DAmP' = SGA_DAmP,
            .id = 'interaction_network') %>%
  separate(`Query Strain ID`, into = c('Query_ORF', 'Query_descriptor'), sep ='_', extra = 'merge', remove = FALSE) %>%
  separate(`Query allele name`, into = c('Query_name', 'Query_mutant'), sep ='-', extra = 'merge', remove = FALSE) %>% 
  separate(`Array Strain ID`, into = c('Array_ORF', 'Array_descriptor'), sep ='_', extra = 'merge', remove = FALSE) %>%
  separate(`Array allele name`, into = c('Array_name', 'Array_mutant'), sep ='-', extra = 'merge', remove = FALSE)

# define function to read and process krogan lab E-MAP screens
load_emap <- function(dataset) {
  
  data <-
    read_delim(paste0('Krogan lab EMAP screens/', dataset, '_gene_names.txt'), delim='\t') %>% 
    rename(Query_strain = Gene) %>% 
    mutate(Query_ID = read_delim(paste0('Krogan lab EMAP screens/', dataset, '.txt'), delim='\t')$Gene) %>%
    select(Query_ID, Query_strain, everything()) %>% 
    separate(Query_ID, into = c('Query_ORF'), sep =' - ', extra = "drop", remove = FALSE) %>%
    separate(Query_strain, into = c('Query_name','Query_perturbation'), sep =' - ', remove = FALSE) %>%
    mutate(Query_perturbation = ifelse(is.na(Query_perturbation), "KO", Query_perturbation))
  
  strain_info <- 
    data %>% 
    select(starts_with('Query')) %>% 
    rename(Array_ID = Query_ID,
           Array_ORF = Query_ORF,
           Array_strain = Query_strain,
           Array_name = Query_name,
           Array_perturbation = Query_perturbation)
  
  data %>%
    gather(-starts_with('Query'), key = 'Array_strain', value = 'score', na.rm = TRUE) %>% 
    left_join(strain_info, by = c('Array_strain' = 'Array_strain')) %>% 
    mutate(dataset = dataset) %>% 
    select(dataset, starts_with('Query'), starts_with('Array'), score)
}

# read in EMAP data
# NOTE: all of these are square matrices, where the rows in "Gene" are the same order as the columns after "Gene"
# additionally, the ORF datafiles will have ORFs in the form "ORF - ORF" for knockdown strains
# (as opposed to "ORF" for knockout strains), and the string is always the same before and after the hyphen.
TF <- load_emap('TF')                     # TF, transcription factor emap [not 1536]
cF3 <- load_emap('cF3')                   # cF3, chromosome biology emap [not 1536]
combRNAP <- load_emap('combRNAP')         # combRNAP, RNA processing [not 1536]
cc2I <- load_emap('cc2I')                 # cc2I, cell cycle screens (unpublished?)
endo <- load_emap('endo')                 # endo, endocytosis emap (unpublished?) [not 1536]
espT <- load_emap('espT')                 # espT, early secretory pathway emap [not 1536]
kinase <- load_emap('kinase')             # kinase, kinase emap [not 1536]
scorematQC1p <- load_emap('scorematQC1p') # scorematQC1p, lipid emap [1536]
slp <- load_emap('slp')                   # slp, mitochondrial emap [1536]

# bind emap data together
EMAP <- bind_rows(TF, cF3, combRNAP, cc2I, endo, espT, kinase, scorematQC1p, slp)

# write datafiles
# Rdata files, well compressed and fast loading
save(SGA, file='SGA_preprocessed.rda')
save(EMAP, file='EMAP_preprocessed.rda')

# Text files, bulky but readable
# write_delim(SGA, 'SGA_preprocessed.txt', delim = '\t')
# write_delim(EMAP, 'EMAP_preprocessed.txt', delim = '\t')

# hdf5 files, lightweight and universal, but slower than Rdata files (at least for these sets)
# library(rhdf5)
# h5write(EMAP, 'EMAP.h5', 'EMAP')
# h5write(SGA, 'SGA.h5', 'SGA')
