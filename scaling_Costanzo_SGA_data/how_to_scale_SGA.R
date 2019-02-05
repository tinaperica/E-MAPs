# getting scaled SGA data
# 
library(tidyverse)
load('data/SGA.rda')

SGA <-
  SGA %>%
  rename('Query_Strain_ID' = `Query Strain ID`, # Clean up to make it easier to type & code
       'Array_Strain_ID' = `Array Strain ID`,
       'sga_score' = `Genetic interaction score (ε)`,
       'pvalue' = `P-value`) %>%
  filter(pvalue < 0.05) %>% 
  #filter(abs(sga_score) > 0.08) %>%                # include for intermediate 
  #filter(sga_score < -0.12 | sga_score > 0.16) %>% # include for stringent
  select(Query_Strain_ID, Query_ORF, Array_Strain_ID, Array_ORF, #ORF1, ORF2,
         sga_score, pvalue)

#load('slp_spline.RData')
load('cF3_spline.RData')

# Scale the
SGA_scaled <-
  SGA %>%
  mutate(scaling_factor = predict(scaling_spline, .$sga_score)$y) %>% 
  mutate(sga_score_scaled = sga_score * scaling_factor)

save(SGA_scaled, file = "cF3_scaled_SGA2016.RData")

                      