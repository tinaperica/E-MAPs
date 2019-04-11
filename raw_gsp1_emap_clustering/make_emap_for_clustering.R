
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

read_delim('gsp1_pEMAP_avg_merged_gene_names.txt', delim = '\t') %>% 
  gather(-Gene, key = strain, value = score) %>% 
  separate(strain, 'library_gene', sep = ' - ', remove = TRUE) %>% 
  separate(Gene, c('Gene', 'mutant'), sep = ' - ') %>% 
  select(mutant, library_gene, score) %>% 
  spread(library_gene, score) %>% 
  arrange(match(mutant, mutants_ordered_by_res)) %>% 
  write_delim('gsp1_emap_for_clustering.txt', delim='\t', na = '')


# Prepare the dataset of core delta rASA in wide format for clustering

core_deltarASA_by_resnum <-
  read_delim('../SASA_interfaces.txt', delim = '\t') %>% 
  rename(resnum = yeastresnum) %>% 
  filter(interface == 'core') %>%
  select(-interface)

core_deltarASA_by_mutant <-
  read_delim('gsp1_mutants_ordered_by_residue.txt', delim="\n", col_names='mutant') %>%
  mutate(resnum = as.numeric(substr(mutant, 2, nchar(mutant)-1))) %>% 
  left_join(core_deltarASA_by_resnum, by = 'resnum') %>% 
  select(-resnum) # %>% ggplot(aes(mutant, partner)) + geom_tile(aes(fill = deltarASA))

core_deltarASA_by_mutant %>% 
  spread(partner, deltarASA, fill=FALSE) %>% 
  select(-`<NA>`) %>% 
  mutate(mutant = factor(mutant, levels = mutants_ordered_by_res)) %>% 
  arrange(mutant) %>% 
  write_delim('core_deltarASA_by_mutant_for_clustering.txt', delim = '\t', na = '')


# Prepare the dataset of GAP/GEF efficiency ratio in wide format (really it's
# just one measurement per mutant) for clustering.

GAP_GEF_ratio <-
  read_delim('../titration_curves/GAP_GEF_ratio.txt', delim = '\t') %>% 
  select(mutant, GAP_GEF_ratio)

GAP_GEF_ratio$mutant[GAP_GEF_ratio$mutant == 'WT'] <- 'GSP1-NAT'

GAP_GEF_ratio %>% 
  mutate(mutant = factor(mutant, levels = mutants_ordered_by_res)) %>% 
  arrange(mutant) %>% 
  write_delim('GAP_GEF_ratio_for_clustering.txt', delim = '\t', na = '')  

