library(pracma)
options(stringsAsFactors = F)
abs.pmax <- function(x, y) {
  df <- data.frame(x, y)
  pmax <- vector()
  for (i in 1:length(x)) {
    pair <- df[i, c("x", "y")]
    max <- pair[ which.max(abs(pair)) ][1,1]
    pmax <- append(pmax, max)
  }
  return(pmax)
}
abs.pmin <- function(x, y) {
  df <- data.frame(x, y)
  pmin <- vector()
  for (i in 1:length(x)) {
    pair <- df[i, c("x", "y")]
    min <- pair[ which.min(abs(pair)) ][1,1]
    pmin <- append(pmin, min)
  }
  return(pmin)
}
pmean <- function(x, y) {
  df <- data.frame(x, y)
  pmean <- vector()
  for (i in 1:length(x)) {
    pair <- df[i, ]
    mean <- mean(pair$x, pair$y, na.rm = T)
    pmean <- append(pmean, mean)
  }
  return(pmean)
}
asym_pnorm_cropped <- function(data) {  ### return 1 - 2*probability that the value is same or less (for neg scores) or same or greater for pos scores
  prob_weights <- vector()
  for (x in data) {
    if (x < -2) {
      if ( x < 0 ) {
        p <- pnorm(x, mean = 0, sd = 5)
     } else {
      p <- 0.5
    }
    } else {
      if (x > 1) {
        p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      } else {
        p <- 0.5
      }
    }
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}
asym_pnorm <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
    } 
    prob_weights <- append(prob_weights, 1 - 2*p)
  }
  return(prob_weights)
}
w_pearson <- function(x, y, w) {
  non_zero_w <- w[w != 0]
  if (length(non_zero_w) > 2) {  #### filter by mean weight or sum of weights - try to empirically get a good threshold
    # x and y weighted means
    mean_x <- sum(x * w) / sum(w)
    mean_y <- sum(y * w) / sum(w)
    # Compute the weighted variance
    vx <- sum( w * (x - mean_x)^2 ) / sum(w)
    vy <- sum( w * (y - mean_y)^2 ) / sum(w)
    # Compute the covariance
    vxy <- sum( w * (x - mean_x) * (y - mean_y) ) / sum(w)
    # Compute the correlation
    w_correlation <- (vxy / sqrt(vx * vy))
    values <- c(w_correlation)
  } else {
    values <- c(NaN)
  }
  return(values)
}
w_pearson_wrap <- function( x, y) {
  ## pick pairwise observable
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y)
  w_cropped <- asym_pnorm_cropped(max)
  w <- asym_pnorm(max)
  df <- cbind(df, data.frame(w, w_cropped))
  df_clean <- df[df$w_cropped != 0, ]
  values <- vector()
  values <- append(values, w_pearson(x, y, w))     # weighted pearson
  values <- append(values, w_pearson(x, y, w_cropped)) # cropped weighted pearson
  if (length(df_clean$x) > 2) {
    values <- append(values, w_pearson(x, y, w))  # weighted pearson, filtered by uncropped points
  } else {
    values <- append(values, NaN)
  }
  return(values)
}
partial_corr <- function(x, y, z) {
  df <- data.frame(x, y, z)
  df <- df[complete.cases(df),]
  lmx <- lm(data = df, x ~ z)
  x.res <- lmx$residuals
  lmy <- lm(data = df, y ~ z)
  y.res <- lmy$residuals
  cor <- cor(x.res, y.res)
  return(cor)
}
#task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
#inputfilename <- paste0("info_for_correlations/", task_n, "_task_info.RData")
#inputfilename <- "clustered_correlations/preprocessed_data_for_correlations/8001066_task_info.RData"
#load(inputfilename)
load("2018-03-09_emap_data_for_corr.RData")
### test weights
sample_x <- sample(ubermap$ubermap$score[! is.na(ubermap$ubermap$score)], 10000)
plot(sample_x, asym_pnorm(sample_x), xlab = "E-MAP score", ylab = "weight")
points(sample_x, asym_pnorm_cropped(sample_x), col = "red")
legend("bottomleft", legend = c("continuous weights", "cropped weights"), fill = c("black", "red"))

#mut_ubermap <- ubermap$ubermap[ubermap$ubermap$ORF == "YLR293C",]
#mean_lib_mut_value <- with(mut_ubermap, aggregate(score, 
 #               by = list(library.gene_name = library.gene_name), na.rm.mean))
#names(mean_lib_mut_value)[2] <- "mean_lib_score" 
#ubermap$ubermap <- merge(ubermap$ubermap, mean_lib_mut_value, by = "library.gene_name")
#outputfilename <- paste0("20180310_", task_n, "_correlations.RData")
#output_file_path <- file.path("Output/correlations_w_random/", outputfilename)
#output_file_path <- file.path("clustered_correlations/correlation_RData/mutants_only.RData") 
correlations_df <- data.frame()
all_weights <- vector()
all_cropped_weights <- vector()
#for ( p in seq_along(task.info[["pairs"]][1,]) ) {
pdf("examples_of_correlations.pdf", width = 10)
for (input in c(1, 8001066)) {
  inputfilename <- paste0("clustered_correlations/preprocessed_data_for_correlations/", input, "_task_info.RData")
  load(inputfilename)
  for(p in c(1:50)) {
    query1 <- task.info[["pairs"]][1, p]
    query2 <- task.info[["pairs"]][2, p]
  
    ubermap.query1 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query1, ]
    ubermap.query2 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query2, ]
   #if (ubermap.query1$ORF[1] != ubermap.query2$ORF[1]) {
    #for ( i in seq_along( ubermap[["clusters"]] ) ) {
    for ( i in 1:10 ) {
  
      cluster <- ubermap[["clusters"]][i]
      temp_library_clusters <- ubermap[["library_clusters"]][ubermap[["library_clusters"]][["cluster"]] == cluster, ]
      temp.ubermap.query1 <- ubermap.query1[ubermap.query1[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    
      if (length(temp.ubermap.query1[,1]) > 2) {   
      
        temp.ubermap.query2 <- ubermap.query2[ubermap.query2[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
      
        if (length(temp.ubermap.query2[,1]) > 2) {
          merged.emap.data <- merge(temp.ubermap.query1[, c("Gene_uniq", "library.ORF","score","random_score", "random_score_2")], 
                              temp.ubermap.query2[, c("Gene_uniq", "library.ORF", "score", "random_score", "random_score_2")], 
                              by = "library.ORF")
          if (length(merged.emap.data[complete.cases(merged.emap.data[, c("score.x", "score.y")]), ][,1]) > 2) {
            correlation <- round(cor(merged.emap.data$score.x, merged.emap.data$score.y, use = "pairwise.complete.obs"), 4)
            wp <- round(w_pearson_wrap(merged.emap.data$score.x, merged.emap.data$score.y), 4)
            w_correlation <- wp[1]
            w_crop_correlation <- wp[2]
            filtered_w_correlation <- wp[3]
            df <- data.frame("x" = merged.emap.data$score.x, 
                             "y" = merged.emap.data$score.y)
            df <- df[complete.cases(df),]
            max <- abs.pmax(df$x, df$y)
            weights <- asym_pnorm(max)
            mean_weight <- round(mean(weights, na.rm = T), 3)
            all_weights <- append(all_weights, weights)
            cropped_weights <- asym_pnorm_cropped(max)
            mean_cropped_weight <- round(mean(cropped_weights, na.rm = T), 3)
            all_cropped_weights <- append(all_cropped_weights, cropped_weights)
            plotting.df <- data.frame("x" = df$x, "y" = df$y, "weight" = cropped_weights)
            fit <- odregress(plotting.df$x, plotting.df$y)
            if (length(plotting.df$x[plotting.df$weight > 0]) > 1) {
              cropped_fit <- odregress(plotting.df$x[plotting.df$weight > 0], plotting.df$y[plotting.df$weight > 0])
            } else {
              cropped_fit$coeff <- c(NaN, NaN)
            }
            slope <- fit$coeff[1]
            cropped_slope <- cropped_fit$coeff[1]
            op <- par(mfrow = c(1, 2))
            plot(plotting.df$x, plotting.df$y, main = cluster, xlab = query1, 
                ylab = query2, xlim = c(-18, 12), ylim = c(-18, 12))
            abline(a = fit$coeff[2], b = fit$coeff[1])
            if (! is.na(cropped_fit$coeff[1])) {
              abline(a = cropped_fit$coeff[2], b = cropped_fit$coeff[1], col = "red")
            }
            points(plotting.df$x[plotting.df$weight > 0], plotting.df$y[plotting.df$weight > 0], col = "red")
            points(plotting.df$x[plotting.df$weight > 0.9], plotting.df$y[plotting.df$weight > 0.9], col = "red", pch = 19)
            legend("top", legend = c(
                    paste0("corr = ", correlation), 
                    paste0("w pearson = " , w_correlation), 
                    paste0("w pearson crop = ", w_crop_correlation),
                    paste0("filtered weighted corr  = ", round(filtered_w_correlation, 4))
                    ),
                    cex = 0.75)
            inverse_fit <- odregress(plotting.df$y, plotting.df$x)
            inverse_slope <- inverse_fit$coeff[1]
            if (length(plotting.df$x[plotting.df$weight > 0]) > 1) {
              cropped_inverse_fit <- odregress(plotting.df$y[plotting.df$weight > 0], plotting.df$x[plotting.df$weight > 0])
            } else {
              cropped_inverse_fit$coeff <- c(NaN, NaN)
            }
            cropped_inverse_slope <- cropped_inverse_fit$coeff[1]
            abs_min_slope <- abs.pmin(slope, inverse_slope)
            abs_max_slope = abs.pmax(slope, inverse_slope)
            abs_min_crop_slope <- abs.pmin(cropped_slope, cropped_inverse_slope)
            abs_max_crop_slope = abs.pmax(cropped_slope, cropped_inverse_slope)
            if (! is.na(w_crop_correlation)) {
              w_crop_corr_min_slope <- w_crop_correlation * abs(abs_min_crop_slope)
              w_crop_corr_max_slope <- w_crop_correlation * abs(abs_max_crop_slope)
            } else {
              w_crop_corr_min_slope <- NaN
              w_crop_corr_max_slope <- NaN
            }
            if (! is.na(filtered_w_correlation)) {
              w_filtered_corr_min_slope <- filtered_w_correlation * abs(abs_min_slope)
              w_filtered_corr_max_slope <- filtered_w_correlation * abs(abs_max_slope)
            } else {
              w_filtered_corr_min_slope <- NaN
              w_filtered_corr_max_slope <- NaN
            }
            corr_min_slope <- correlation * abs(abs_min_slope)
            corr_max_slope <- correlation * abs(abs_max_slope)
            plot(plotting.df$y, plotting.df$x, xlab = query2,
            ylab = query1, xlim = c(-18, 12), ylim = c(-18, 12))
            abline(a = inverse_fit$coeff[2], b = inverse_fit$coeff[1])
            if (! is.na(cropped_inverse_slope)) {
            abline(a = cropped_inverse_fit$coeff[2], b = cropped_inverse_slope, col = "red")
            }
            points(plotting.df$y[plotting.df$weight > 0], plotting.df$x[plotting.df$weight > 0], col = "red")
            points(plotting.df$y[plotting.df$weight > 0.9], plotting.df$x[plotting.df$weight > 0.9], col = "red", pch = 19)
            legend("top", legend = c(
                paste0("slope = ", round(abs_min_crop_slope, 2)),
                paste0("corr * min slope = ", round(corr_min_slope, 2)),
                paste0("weighted crop corr * min slope = ", round(w_crop_corr_min_slope, 3)),
                paste0("filtered weighted corr * min slope = ", round(w_filtered_corr_min_slope, 2))
            ),
            cex = 0.75)
          #} else {
            #correlation <- NaN
            #w_correlation <- NaN
            #spearman <- NaN
            #w_crop_correlation <- NaN
            #partial_correlation <- NaN
          #}
          # if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score.x", "random_score.y")]), ][,1]) > 5) {
          #   random_correlation <- round(cor(merged.emap.data$random_score.x, merged.emap.data$random_score.y, use =  "pairwise.complete.obs"), 5)
          #   w_random_correlation <- round(w_pearson(merged.emap.data$random_score.x, merged.emap.data$random_score.y), 5)
          # } else {
          #   random_correlation <- "NA"
          #   w_random_correlation <- "NA"
          # }
          # # if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score_2.x", "random_score_2.y")]), ][,1]) > 5) {
          #   random_correlation_2 <- round(cor(merged.emap.data$random_score_2.x, merged.emap.data$random_score_2.y, use =  "pairwise.complete.obs"), 5)
          #   w_random_correlation_2 <- round(w_pearson(merged.emap.data$random_score_2.x, merged.emap.data$random_score_2.y), 5)
          # } else {
          #   random_correlation_2 <- "NA"
          #   w_random_correlation_2 <- "NA"
          # }
          # if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score.x", "random_score_2.y")]), ][,1]) > 10) {
          #   random_correlation_3 <- round(cor(merged.emap.data$random_score.x, merged.emap.data$random_score_2.y, use = "pairwise.complete.obs"), 3)
          # } else {
          #   random_correlation_3 <- "NA"
          # }
          # if ( length(merged.emap.data[complete.cases(merged.emap.data[, c("random_score_2.x", "random_score.y")]), ][,1]) > 10) {
          #   random_correlation_4 <- round(cor(merged.emap.data$random_score_2.x, merged.emap.data$random_score.y, use = "pairwise.complete.obs"), 3)
          # } else {
          #   random_correlation_4 <- "NA"
          }
          
          correlations_df <- rbind(correlations_df, data.frame( 
                "Gene_uniq1" = query1, "Gene_uniq2" = query2, "cluster" = cluster, 
                correlation, w_correlation, w_crop_correlation, filtered_w_correlation, "slope1" = slope, "slope2" = inverse_slope, 
                "max_slope" = max(slope, inverse_slope, na.rm = T), "min slope" = min(slope, inverse_slope, na.rm = T),
                abs_max_slope, abs_min_slope,
                w_crop_corr_min_slope, w_crop_corr_max_slope, 
                w_filtered_corr_min_slope, w_filtered_corr_max_slope, 
                corr_min_slope, corr_max_slope,
                mean_weight, mean_cropped_weight
                ))
        }
      }
    }
   }
  #}
}
dev.off()
op <- par(mfrow = c(1, 1))

head(correlations_df)

correlations_df <- correlations_df[! is.na(correlations_df$Gene_uniq1), ]
head(correlations_df[! is.na(correlations_df$w_crop_correlation) & is.na(correlations_df$filtered_w_correlation),])

##### these plots are based on example 1e5 calculations
## saved in example_corr_df.RData  
correlations_df <- correlations_df[order(correlations_df$w_crop_corr_min_slope),]
plot(density(correlations_df$corr_min_slope, na.rm = T), col = "orange", lwd = 2, main = "Distribution of values")
lines(density(correlations_df$w_crop_corr_min_slope, na.rm = T), col = "blue", main = "Distribution of correlation measures")
lines(density(correlations_df$correlation), lwd = 2)
lines(density(correlations_df$w_correlation), col = "green", lwd = 2)
lines(density(correlations_df$w_crop_correlation, na.rm = T), col = "red", lwd = 2)
legend("topleft", legend = c("corr min slope", "w crop corr min slope", "corr", "w corr", "w crop corr"),
       fill = c("orange", "blue", "black", "green", "red"))
plot(density(correlations_df$correlation, na.rm = T), lwd = 2, main = "Distribution of values")
lines(density(correlations_df$w_correlation, na.rm = T), col = "green", lwd = 2)
lines(density(correlations_df$w_crop_correlation, na.rm = T), col = "red", lwd = 2)
legend("topleft", legend = c("Pearson correlation", "weighted Pearson", "weighted cropped Pearson"),
       fill = c("black", "green", "red"))
hist(correlations_df$w_crop_corr_min_slope)


cor <- round(cor(correlations_df$w_correlation, correlations_df$correlation), 3)
fit <- lm(data = correlations_df, w_correlation ~ correlation)
plot(correlations_df$correlation, correlations_df$w_correlation,
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05), pch = 20,
     xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Pearson versus weighted Pearson", 
     xlab = "Pearson corr", ylab = "Pearson weighted corr")
abline(b = 1, a = 0)
abline(fit, col = "blue", lwd = 2)
legend("topleft", legend = c(paste("corr = ", cor), 
        paste("slope = ", round(fit$coefficients[2], 2))))

cor <- round(cor(correlations_df$w_crop_correlation, correlations_df$correlation, use = "pairwise.complete.obs"), 3)
fit <- lm(data = correlations_df, w_crop_correlation ~ correlation)
plot(correlations_df$correlation, correlations_df$w_crop_correlation, 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.02), pch = 20,
     xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Pearson versus cropped weighted Pearson", 
     xlab = "Pearson corr", ylab = "Pearson cropped weighted corr")
abline(b = 1, a = 0)
abline(fit, col = "blue", lwd = 2)
legend("topleft", legend = c(paste("corr = ", cor), 
        paste("slope = ", round(fit$coefficients[2], 2))))


cor <- round(cor(correlations_df$w_crop_corr_min_slope, correlations_df$correlation, use = "pairwise.complete.obs"), 3)
fit <- lm(data = correlations_df,  w_crop_corr_min_slope ~ correlation)
plot(correlations_df$correlation, correlations_df$w_crop_corr_min_slope, 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.02), pch = 20,
     xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Pearson versus cropped weighted slope scaled Pearson", 
     xlab = "Pearson corr", ylab = "Pearson cropped weighted slope scaled corr")
abline(b = 1, a = 0)
abline(fit, col = "blue", lwd = 2)
legend("topleft", legend = c(paste("corr = ", cor), 
                paste("slope = ", round(fit$coefficients[2], 2))))

cor <- round(cor(correlations_df$w_filtered_corr_min_slope, correlations_df$correlation, use = "pairwise.complete.obs"), 3)
fit <- lm(data = correlations_df, w_filtered_corr_min_slope ~ correlation)
plot(correlations_df$correlation, correlations_df$w_filtered_corr_min_slope,
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.02), pch = 20, 
     xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Pearson versus filtered weighted slope scaled Pearson", 
     xlab = "Pearson corr", ylab = "Pearson filtered weighted slope scaled corr")
abline(b = 1, a = 0)
abline(fit, col = "blue", lwd = 2)
legend("topleft", legend = c(paste("corr = ", cor), 
                paste("slope = ", round(fit$coefficients[2], 2))))


cor <- round(cor(correlations_df$filtered_w_correlation, correlations_df$correlation, use = "pairwise.complete.obs"), 3)
fit <- lm(data = correlations_df, filtered_w_correlation ~ correlation)
plot(correlations_df$correlation, correlations_df$filtered_w_correlation,
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.02), pch = 20, 
     xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Pearson versus filtered weighted Pearson", 
     xlab = "Pearson corr", ylab = "Pearson filtered weighted corr")
abline(b = 1, a = 0)
abline(fit, col = "blue", lwd = 2)
legend("topleft", legend = c(paste("corr = ", cor), paste("slope = ", round(fit$coefficients[2], 2))))


head(correlations_df[correlations_df$w_crop_corr_min_slope > 0.5,])
#scaling.fit <- lm(data = correlations_df, scaled_spread_weight ~ spread_weight)
# Call:
#   lm(formula = scaled_spread_weight ~ spread_weight, data = correlations_df)
# 
# Coefficients:
#   (Intercept)  spread_weight  
# -0.76894        0.04122  
# 
#plot(correlations_df$correlation, correlations_df$w_correlation)
#abline(a = 0, b = 1 )
#plot(correlations_df$correlation, correlations_df$w_crop_correlation)
#abline(a = 0, b = 1 )
#correlations_df <- cbind(correlations_df, "w_ratio" = correlations_df$w_crop_correlation/correlations_df$w_correlation)
#head(correlations_df[correlations_df$w_crop_correlation > 0.9, ])
plot(correlations_df$w_correlation, correlations_df$w_crop_correlation)
abline(a = 0, b = 1)

plot(density(as.numeric(correlations_df$correlation), na.rm = T))
lines(density(as.numeric(correlations_df$w_correlation), na.rm = T), col = "red")
lines(density(as.numeric(correlations_df$partial_correlation), na.rm = T), col = "green")
lines(density(as.numeric(correlations_df$w_crop_correlation), na.rm = T), col = "blue")
plot(density(weights))
na.rm.mean <- function(x) {
  mean(x, na.rm = T)
} 

### try partial correlation 
# for Z use mean value score for each library gene (try mean for the whole ubermap (might be zero)
# or mean for mutants)
#mean_lib_value <- with(ubermap$ubermap, aggregate(score, 
              #by = list(library.gene_name = library.gene_name), na.rm.mean))                       
#mean_lib_mut_value <- with(mut_ubermap, aggregate(score, 
 #           by = list(library.gene_name = library.gene_name), na.rm.mean))


w_vs_crop_w <- correlations_df[(! is.na(correlations_df$w_correlation)) & is.na(correlations_df$w_crop_correlation),]
#plot(all_weights, all_cropped_weights)
plot(density(correlations_df$mean_weight))
lines(density(correlations_df$mean_cropped_weight), col = "blue")

plot(correlations_df$w_crop_corr_min_slope, correlations_df$w_crop_corr_max_slope)
plot(correlations_df$corr_min_slope, correlations_df$corr_max_slope)
