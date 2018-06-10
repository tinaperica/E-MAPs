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
  } else {
    w_correlation <- NaN
  }
  return(w_correlation)
}
correlation_wrap <- function( x, y) {
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y)
  w <- asym_pnorm_cropped(max)
  df <- cbind(df, w)
  corr <- w_pearson(x, y, w)
  if (length(df$x[df$w > 0]) > 1) {
    fit <- odregress(df$x[df$w > 0], df$y[df$w > 0])
    inverse_fit <- odregress(df$y[df$w > 0], df$x[df$w > 0])
  } else {
    fit$coeff <- c(NaN, NaN)
    inverse_fit$coeff <- c(NaN, NaN)
  }
  slope <- fit$coeff[1]
  inverse_slope <- inverse_fit$coeff[1]
  
  
}

task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
inputfilename <- paste0("info_for_correlations/", task_n, "_task_info.RData")
load(inputfilename)
load("2018-03-09_emap_data_for_corr.RData")

for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  ubermap.query1 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query1, ]
  ubermap.query2 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query2, ]
  for ( i in seq_along( ubermap[["clusters"]] ) ) {
    cluster <- ubermap[["clusters"]][i]
    temp_library_clusters <- ubermap[["library_clusters"]][ubermap[["library_clusters"]][["cluster"]] == cluster, ]
    temp.ubermap.query1 <- ubermap.query1[ubermap.query1[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    temp.ubermap.query2 <- ubermap.query2[ubermap.query2[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    merged.emap.data <- merge(temp.ubermap.query1[, c("Gene_uniq", "library.ORF","score","random_score", "random_score_2")], 
                              temp.ubermap.query2[, c("Gene_uniq", "library.ORF", "score", "random_score", "random_score_2")], 
                              by = "library.ORF")
    if (length(merged.emap.data[complete.cases(merged.emap.data[, c("score.x", "score.y")]), ][,1]) > 2) {
        w_crop_correlation <- round(w_pearson_wrap(merged.emap.data$score.x, merged.emap.data$score.y), 4)
        df <- data.frame("x" = merged.emap.data$score.x, 
                           "y" = merged.emap.data$score.y)
        df <- df[complete.cases(df),]
          max <- abs.pmax(df$x, df$y)
          plotting.df <- data.frame("x" = df$x, "y" = df$y, "weight" = cropped_weights)
          if (length(plotting.df$x[plotting.df$weight > 0]) > 1) {
            cropped_fit <- odregress(plotting.df$x[plotting.df$weight > 0], plotting.df$y[plotting.df$weight > 0])
          } else {
            cropped_fit$coeff <- c(NaN, NaN)
          }
          cropped_slope <- cropped_fit$coeff[1]
          if (length(plotting.df$x[plotting.df$weight > 0]) > 1) {
            cropped_inverse_fit <- odregress(plotting.df$y[plotting.df$weight > 0], plotting.df$x[plotting.df$weight > 0])
          } else {
            cropped_inverse_fit$coeff <- c(NaN, NaN)
          }
          cropped_inverse_slope <- cropped_inverse_fit$coeff[1]
          
          abs_min_crop_slope <- abs.pmin(cropped_slope, cropped_inverse_slope)
          if (! is.na(w_crop_correlation)) {
            w_crop_corr_min_slope <- w_crop_correlation * abs(abs_min_crop_slope)
          } else {
            w_crop_corr_min_slope <- NaN
          }
        }
      }
    }
        
      correlations_df <- rbind(correlations_df, data.frame( 
          "Gene_uniq1" = query1, "Gene_uniq2" = query2, "cluster" = cluster, 
          correlation, w_correlation, w_crop_correlation, "slope1" = slope, "slope2" = inverse_slope, 
          "max_slope" = max(slope, inverse_slope, na.rm = T), "min slope" = min(slope, inverse_slope, na.rm = T),
          abs_max_slope, abs_min_slope,
          w_crop_corr_min_slope, w_crop_corr_max_slope, corr_min_slope, corr_max_slope,
          mean_weight, mean_cropped_weight
        ))
    }
  }
}


