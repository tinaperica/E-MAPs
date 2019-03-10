# NOTE This script assumes any averaging for EMAP data and SGA data has already been done.

##### Pseudo-code for SGA preprocessing #####
#  - Pre-process EMAP and SGA data: choose reference EMAP dataset, do any averaging you want
#     - "First, in order to have a single measurement for each SGA interaction, we average
#       all interactions where multiple measurements were available (eg query geneA x array 
#       geneB and query geneB array geneA). For interactions with only a single measurement
#       available, a pseudo-averaging was performed, following...Collins et al 2006"
#     - Maybe don't need to do this, write the code for it anyway, and try both.
#  - Scale the average SGA scores to make them comparable with the reference E-MAP data:
#     - Take pairs interactions in both SGA and EMAP reference
#     - Order each subset from lowest to highest
#     - Split into 500 bins for each subset
#         - 2nd bin has 50% overlap with first and third, and so on
#     - Compare mean of each SGA bin to mean of EMAP bin, computing scaling factor as the
#       ratio mu_EMAP/mu_SGA
#     - Fit a univariate spline to the function scaling factor = f(SGA score)
#         - Do not use SGA scores close to zero to fit this curve
#     - A scaling factor for each SGA score in the entire data subset can be caluclated using
#       this fitted curve
#   - Compare the QQ plots from before and after scaling. In Ryan et al 2012, the scaled
#     SGA data was nearly identical to the EMAP data in the range of -10 to 5. Points outside
#     of this range correspond to less than 0.1% of the data.
#   - This method requires assuming that (1) there exists an identical number of "true"
#     genetic interactions shared between the two datasets, and (2) that the rank ordering of
#     both methods is of similar quality. Ryan et al 2012 asserts (1) is trivially true and
#     showed (2) by using EMAP and SGA scores to predict biologically known synthetic lethal
#     interactions from BIOGRID.
#   - Comparison of EMAP datasets suggests that they're highly correlated but not identical.

##### Import Packages #####
library(tidyverse)

##### Load Data #####
load('data/SGA.rda')
load('data/EMAP.rda')
load('data/June2016_Gsp1_E-MAP_data.RData')

##### Choose a Reference EMAP Dataset to use ##### 
# Based on the largest overlap with our Gsp1 EMAP (it's the slp dataset)
# FIX -- this doesnt answer the question I say it does???
# select(EMAP, Array_ORF, dataset) %>% unique() %>% count(by=dataset)
# count(e.map, by=mutant)

##### Prepare the EMAP Reference Dataset #####
# We need to remove the duplicate measurements for the slp dataset The slp dataset was square, with
#   repeat measurements for each value, so that the score for pair (Q,A) is also recorded for pair (A,Q).
# This can be shown by ordering the perturbations in a pair in a standard way (alphabetically) and then
#   grouping by the standard ordering (Pert1, Pert2). To confirm each group contains exactly two
#   duplicates, call:
#     group_by(Pert1, Pert2) %>% mutate(n=n(), sd = sd(score)) %>% summary()
# There are 869,676 observations in the slp dataset, and 434,838 observations after removing duplicates.
#   To confirm that no two double perturbation strains are repeated in this cleaned dataset, call:
#   group_by(Query_strain, Array_strain) %>% filter(n() > 1)

EMAP_ref <-
  EMAP %>%
  filter(dataset == 'cF3') %>% 
  rowwise() %>% 
  mutate(Pert1 = min(Query_strain, Array_strain),
         Pert2 = max(Query_strain, Array_strain)) %>%
  ungroup() %>% 
  distinct(Pert1, Pert2, .keep_all = TRUE) %>% 
  rename(emap_score = score) %>% 
  select(-starts_with('Pert'))

##### Prepare the SGA Dataset #####
# Both the EMAP and SGA datasets contain multiple distinct strains targeting the same gene. At
#   this stage we want to keep all distinct strains, so we'll intersect the datasets based on
#   the unique ORFs
EMAP_ORF_pairs <- 
  EMAP_ref %>% 
  select(ends_with('ORF')) %>% 
  rowwise() %>% 
  mutate(ORF1 = min(Query_ORF, Array_ORF),
         ORF2 = max(Query_ORF, Array_ORF)) %>%
  select(starts_with('ORF')) %>% 
  ungroup() %>% 
  distinct()

# Costanzo et al (2016) recommends three confidence thresholds for using SGA scores: lenient
#   (p<0.05), intermediate(p<0.05 & |eps|>0.08), and stringent (p<0.05 and eps>0.16 | eps<-0.12)
# RUNTIME: takes about 13min

SGA_overlap <-
  SGA %>%
  rename('Query_Strain_ID' = `Query Strain ID`, # Clean up to make it easier to type & code
       'Array_Strain_ID' = `Array Strain ID`,
       'sga_score' = `Genetic interaction score (ε)`,
       'arraytype_temp' = `Arraytype/Temp`,
       'pvalue' = `P-value`) %>%
  filter(pvalue < 0.05) %>% 
  # filter(abs(sga_score) > 0.08) %>%                # include for intermediate 
  # filter(sga_score < -0.12 | sga_score > 0.16) %>% # include for stringent
  rowwise() %>%
  mutate(ORF1 = min(Query_ORF, Array_ORF),
         ORF2 = max(Query_ORF, Array_ORF)) %>%
  ungroup() %>%
  semi_join(EMAP_ORF_pairs, by = c('ORF1', 'ORF2')) %>%
  select(Query_Strain_ID, Query_ORF, Array_Strain_ID, Array_ORF, ORF1, ORF2,
         sga_score, arraytype_temp, pvalue)


# Duplicate measurements: There are duplicate measurements for some double perturbations at two
#  temperatures (DMA30 and TSA26). If we only consider interactions with pvalue < 0.05, this
#  constitutes 216 observations (108 pairs) out of the 79,690 significant SGA observations that
#  overlap with our reference slp EMAP dataset (0.27%)mand 58 obersvations of the 40527
#  significant SGA observations that overlap with our reference cF3 dataset (.14%). These are
#  not enough to warrant concern for skewing the scaling, especially given that these are
#  biologically distinct replicates, just as two different mutants at the same temperature would be.
SGA_overlap %>%
  filter(pvalue < 0.05) %>%
  select(Query_Strain_ID, Array_Strain_ID, arraytype_temp, sga_score) %>%
  rowwise() %>%
  mutate(Pert1 = min(Query_Strain_ID, Array_Strain_ID),
         Pert2 = max(Query_Strain_ID, Array_Strain_ID)) %>%
  group_by(Pert1, Pert2) %>%
  mutate(n_rep = n(), sd_rep = sd(sga_score), mean_rep = mean(sga_score)) %>%
  ungroup() %>%
  filter(n_rep > 1) %>%
  arrange(Pert1, Pert2)

##### Perform the scaling #####

# Get an equal number of scores from each distribution and order them

set.seed(1)
SGA_scores <-
  SGA_overlap %>% 
  sample_n(min(nrow(EMAP_ref),nrow(SGA_overlap)), replace=FALSE) %>% 
  select(sga_score) %>% arrange(sga_score)
EMAP_scores <-
  EMAP_ref %>% 
  sample_n(min(nrow(EMAP_ref),nrow(SGA_overlap)), replace=FALSE) %>% 
  select(emap_score) %>% arrange(emap_score)

# Prepare constants for binning the data. For bins to overlap by 50%, we need 2n-1
# overlapping bins to cover the same subset asn non-overlapping bins.
n_scores <- length(EMAP_scores$emap_score)
n_bins <- 500
bin_interval = floor(n_scores/(n_bins-1))

# bin the scores, then compute means for each bin, and a scaling factor (ratio of means)
binned_data <- data.frame(bin = seq(1:n_bins),
                          start = bin_interval * seq(1:n_bins) - bin_interval,
                          end = bin_interval * seq(1:n_bins) + bin_interval)
binned_data$end[binned_data$bin == 500] <- n_scores # edge case for not-perfectly-divisible n_scores

binned_data <-
  binned_data %>% 
  rowwise() %>% 
  mutate(SGA_mean = mean(SGA_scores$sga_score[start:end], na.rm=T),
         EMAP_mean = mean(EMAP_scores$emap_score[start:end], na.rm=T)) %>% 
  mutate(scaling_factor = EMAP_mean/SGA_mean)
ggplot(binned_data, aes(x=EMAP_mean,y=SGA_mean)) + geom_point()
ggplot(binned_data, aes(x=SGA_mean,y=scaling_factor)) + geom_point()
ggplot(binned_data, aes(x=SGA_mean,y=scaling_factor)) + geom_point() + xlim(-0.1,0.1) + ylim(-10, 10)


# fix binning for slp dataset
# we can see that the values close to zero diverge because dividing by the means
# gives unreliable results, as described in Ryan et al (2012). Accordingly, we find the
# closest critical points near zero and only fit bins that lie outside of those thresholds
lower_sf <- min(filter(binned_data, SGA_mean < 0)$scaling_factor)
lower_bin <- filter(binned_data, scaling_factor == lower_sf)$bin
upper_bin <- filter(binned_data, SGA_mean < 0.07 & SGA_mean > 0, abs(scaling_factor - lower_sf) < 0.005)$bin

scaling_spline <- with(filter(binned_data, bin < lower_bin | bin > upper_bin),
                       smooth.spline(SGA_mean, scaling_factor))
save(scaling_spline, file="slp_spline.RData")
df_spline <- data.frame(predict(scaling_spline, seq(-0.75, 0.5, 0.01)))
ggplot(binned_data, aes(x = SGA_mean, y = scaling_factor)) + geom_point() +
  geom_line(data = df_spline, mapping = aes(x=x, y=y, color='red'))

# fix binning for cF3 dataset
# we can see that the values close to zero diverge because dividing by the means
# gives unreliable results, as described in Ryan et al (2012). Accordingly, we find the
# closest critical points near zero and only fit bins that lie outside of those thresholds
upper_sf <- min(filter(binned_data, SGA_mean > 0)$scaling_factor)
upper_bin <- filter(binned_data, scaling_factor == upper_sf)$bin
lower_bin <- filter(binned_data, SGA_mean < 0, abs(scaling_factor - upper_sf) < 0.005)$bin

scaling_spline <- with(filter(binned_data, bin < lower_bin | bin > upper_bin),
                       smooth.spline(SGA_mean, scaling_factor))
save(scaling_spline, file="cF3_spline.RData")
df_spline <- data.frame(predict(scaling_spline, seq(-0.75, 0.5, 0.01)))
ggplot(binned_data, aes(x = SGA_mean, y = scaling_factor)) + geom_point() +
  geom_line(data = df_spline, mapping = aes(x=x, y=y, color='red'))

# Prep a clean EMAP set
EMAP_clean <- 
  EMAP_ref %>% 
  select(ends_with('ORF'), emap_score) %>% 
  rowwise() %>% 
  mutate(ORF1 = min(Query_ORF, Array_ORF),
         ORF2 = max(Query_ORF, Array_ORF)) %>% 
  ungroup() %>% 
  select(ORF1, ORF2, emap_score)

# Scale the SGA_overlap
SGA_scaled <-
  SGA_overlap %>%
  select(starts_with('ORF'), sga_score) %>% 
  mutate(scaling_factor = predict(scaling_spline, .$sga_score)$y) %>% 
  mutate(sga_score_scaled = sga_score * scaling_factor)

# also try simply scaling to max/min of emap
max_emap <- max(EMAP_clean$emap_score); min_emap <- min(EMAP_clean$emap_score)
max_sga <- max(SGA_scaled$sga_score); min_sga <- min(SGA_scaled$sga_score)

SGA_scaled <-
  SGA_scaled %>% 
  mutate(sga_score_maxmin_scaled = (sga_score - min_sga)/(max_sga - min_sga)*(max_emap - min_emap) + min_emap)

# also try scaling to the regression through the unscaled qqplot
qq_unscaled <- qqplot(x=SGA_EMAP_join$sga_score, y=SGA_EMAP_join$emap_score)
model <- lm(qq_unscaled$y ~ qq_unscaled$x)
coeffs <- model$coefficients

SGA_scaled <-
  SGA_scaled %>% 
  mutate(sga_score_qq_scaled = sga_score*coeffs[2]+coeffs[1])

# plot SGA vs EMAP for overlapping points before and after scaling
SGA_EMAP_join <- inner_join(SGA_scaled, EMAP_clean)

##### Plots #####

# plot of unscaled SGA histogram  
ggplot(SGA_EMAP_join, aes(x=sga_score)) + geom_histogram(binwidth = 0.01)
ggsave('cF3_output/cF3_sga_unscaled_hist.png')

# plot of scaled SGA histogram
ggplot(SGA_EMAP_join, aes(x=sga_score_scaled)) + geom_histogram(binwidth = 0.01)
ggsave('cF3_output/cF3_sga_scaled_hist.png')

# plot of slp EMAP histogram
ggplot(SGA_EMAP_join, aes(x=emap_score)) + geom_histogram(binwidth = 0.01)
ggsave('cF3_output/cF3_emap_hist.png')

# plot of unscaled SGA vs EMAP
ggplot(SGA_EMAP_join, aes(x = emap_score, y = sga_score)) + geom_point()
ggsave('cF3_output/cF3_sga_emap_unscaled.png')

# plot of scaled SGA vs EMAP
ggplot(SGA_EMAP_join, aes(x = emap_score, y = sga_score_scaled)) + geom_point()
ggsave('cF3_output/cF3_sga_emap_scaled.png')

# plot of max/min scaled SGA vs EMAP
ggplot(SGA_EMAP_join, aes(x = emap_score, y = sga_score_maxmin_scaled)) + geom_point()
ggsave('cF3_output/cF3_sga_emap_minmax_scaled.png')

# plot of qq scaled SGA vs EMAP
ggplot(SGA_EMAP_join, aes(x = emap_score, y = sga_score_qq_scaled)) + geom_point()
ggsave('cF3_output/cF3_sga_emap_qq_scaled.png')

# unscaled qqplot
png('cF3_output/cF3_qq_unscaled.png')
qq_unscaled <- qqplot(x=SGA_EMAP_join$emap_score, y=SGA_EMAP_join$sga_score)
model <- lm(qq_unscaled$y ~ qq_unscaled$x)
model$coefficients
plot(qq_unscaled$x, qq_unscaled$y); abline(model$coefficients)
dev.off()

# scaled qqplot
png('cF3_output/cF3_qq_scaled.png')
z <- qqplot(x=SGA_EMAP_join$emap_score, y=SGA_EMAP_join$sga_score_scaled); abline(c(0,1))
dev.off()

# max/min scaled qqplot
png('cF3_output/cF3_qq_maxmin_scaled.png')
z <- qqplot(x=SGA_EMAP_join$emap_score, y=SGA_EMAP_join$sga_score_maxmin_scaled); abline(c(0,1))
dev.off()

# qq scaled qqplot
png('cF3_output/cF3_qq_qq_scaled.png')
z <- qqplot(x=SGA_EMAP_join$emap_score, y=SGA_EMAP_join$sga_score_qq_scaled); abline(c(0,1))
dev.off()
