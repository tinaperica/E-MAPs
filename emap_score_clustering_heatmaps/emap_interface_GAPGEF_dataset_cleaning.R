
library(tidyverse)

# For consistency in visualization of dendrograms, always order mutants by
# residue number, and by amino acid group for mutants of the same residue.
# Groups used in (arbitrary) order:
#   Hydrophobic: AVILMFYW
#   Special cases (we only have G): GPCU
#   Polar uncharged: STNQ
#   Negatively charged: DE
#   Positively charged: RHK

mutants_ordered_by_res <- scan('gsp1_mutants_ordered_by_residue.txt',
                               what="", sep="\n")

# Re-organize the Gsp1 EMAP so that mutants are no longer in form 'GSP1 - A180T'
# but rather 'A180T', and so that the array gene DAmPs are no longer in the form
# 'TAF4 - DAmP' but rather 'TAF4'. Save this as a new txt file for use

read_delim('datasets/gsp1_pEMAP_avg_merged_gene_names.txt', delim = '\t') %>% 
  gather(-Gene, key = strain, value = score) %>% 
  separate(strain, 'library_gene', sep = ' - ', remove = TRUE) %>% 
  separate(Gene, c('Gene', 'mutant'), sep = ' - ') %>% 
  select(mutant, library_gene, score) %>% 
  spread(library_gene, score) %>% 
  arrange(match(mutant, mutants_ordered_by_res)) %>% 
  write_delim('datasets/gsp1_emap_for_clustering.txt', delim='\t', na = '')


# Prepare the dataset of core delta rASA in wide format for clustering

core_deltarASA_by_resnum <-
  read_delim('../../ran_structures/SASA_interfaces.txt', delim = '\t') %>% 
  rename(resnum = yeastresnum) %>% 
  select(-interface)

core_deltarASA_by_mutant <-
  read_delim('gsp1_mutants_ordered_by_residue.txt', delim="\n", col_names='mutant') %>%
  mutate(resnum = as.numeric(substr(mutant, 2, nchar(mutant)-1))) %>% 
  left_join(core_deltarASA_by_resnum, by = 'resnum') %>% 
  drop_na() %>%
  select(-resnum) %>%
  spread(partner, deltarASA, fill=FALSE) %>%
  mutate(mutant = factor(mutant, levels = mutants_ordered_by_res)) %>%
  arrange(mutant) %>% 
  write_delim('datasets/core_deltarASA_by_mutant_for_clustering.txt', delim = '\t')

# Prepare the dataset of GAP and GEF parameters

GAP_GEF_ratio <-
  read_delim('../titration_curves/GAP_GEF_ratio.txt', delim = '\t') %>%
  mutate('ln_GAP_kcat_Km' = log(GAP_kcat_Km),
         'ln_GEF_kcat_Km' = log(GEF_kcat_Km))

GAP_GEF_ratio$mutant[GAP_GEF_ratio$mutant == 'WT'] <- 'GSP1-NAT'

GAP_GEF_ratio %>% 
  mutate(mutant = factor(mutant, levels = mutants_ordered_by_res)) %>% 
  arrange(mutant) %>% 
  write_delim('datasets/GAP_GEF_for_clustering.txt', delim = '\t', na = '') 


# Prepare the Mass Spec ordering
read_delim('apms_log2_fold_change.txt', delim = '\t', col_types = cols()) %>% 
  filter(gene_name %in% c('SRM1', 'RNA1'), norm =='eqM') %>% 
  select(mutant, sample, gene_name, log2FC) %>% 
  spread(gene_name, log2FC) %>% 
  mutate('delta_GAP_GEF' = RNA1 - SRM1) %>% 
  group_by(mutant) %>% 
  mutate('delta_GAP_GEF' = mean(delta_GAP_GEF)) %>%
  ungroup() %>%
  select(mutant, delta_GAP_GEF) %>%
  unique() %>%
  arrange(delta_GAP_GEF) %>%
  select(mutant) %>%
  write_delim('ms_ordering.txt', col_names = FALSE)
