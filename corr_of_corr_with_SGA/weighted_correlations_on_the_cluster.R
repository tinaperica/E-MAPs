remove(list = ls())
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
asym_pnorm_signed <- function(data) {  
  prob_weights <- vector()
  for (x in data) {
    if (x < 0) {
      p <- pnorm(x, mean = 0, sd = 5)
      prob_weights <- append(prob_weights, -1 * (1 - 2*p))
    } else {
      p <- pnorm(x, mean = 0, sd = 2.5, lower.tail = F)
      prob_weights <- append(prob_weights, 1 - 2*p)
    } 
  }
  return(prob_weights)
}
w_pearson <- function( x, y) {
  ## pick pairwise observable
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  x <- df$x
  y <- df$y
  max <- abs.pmax(x, y) # max absolute value E-MAP score of the two
  w_cropped <- asym_pnorm_cropped(max)  ### 0 for E-MAP scores around 0, used for filtering
  w <- asym_pnorm(max)   ## get the weight for the pair
  df <- cbind(df, data.frame(w, w_cropped))
  df_clean <- df[df$w_cropped != 0, ]  ## df_clean is only for filtering
  if (length(df_clean$x) > 1 & length(df$x) > 3) {
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
SoftCosSim <- function(x, y) {
  df <- data.frame(x, y)
  df <- df[complete.cases(df),]
  df$wx <- asym_pnorm_signed(df$x)
  df$wy <- asym_pnorm_signed(df$y)
  max <- abs.pmax(df$x, df$y)
  w_cropped <- asym_pnorm_cropped(max)
  df <- cbind(df, data.frame(w_cropped))
  df_clean <- df[df$w_cropped != 0, ]
  df$w_sim <- (2 - abs(df$wx - df$wy)) / 2
  x <- df$x
  y <- df$y
  w_sim <- df$w_sim
  if (length(df_clean$x) > 1 & length(df$x) > 3) {
    soft_cos_sim <-   sum(w_sim * x * y) / ( sqrt( sum(w_sim * x^2) ) * sqrt( sum(w_sim * y^2) ) )
  } else {
    soft_cos_sim <- NaN
  }
  return(soft_cos_sim)
}
compare_wrap <- function(x, y) {
  values <- vector()
  values <- append(values, w_pearson(x, y))
  values <- append(values, SoftCosSim(x, y))
  return(values)
}

task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
inputfilename <- paste0("all_correlations_task_info/", task_n, "_task_info.RData")
outputfilename <- paste0("20190101_", task_n, "_correlations.RData")
load(inputfilename)
load("20181230_ubermap_for_correlations.RData")
output_file_path <- file.path("Output/correlations", outputfilename)
correlations_df <- data.frame()
for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  ubermap.query1 <- ubermap[ubermap[["query_uniq"]] == query1, ]
  ubermap.query2 <- ubermap[ubermap[["query_uniq"]] == query2, ]
  merged.emap.data <- merge(ubermap.query1, ubermap.query2, by = "library_ORF")
  merged.emap.data <- merged.emap.data[complete.cases(merged.emap.data),]
  values <- round(compare_wrap(merged.emap.data$score.x, merged.emap.data$score.y), 4)
  correlations_df <- rbind(correlations_df, data.frame(
    "query_uniq1" = query1, "query_uniq2" = query2,
    "w_correlation" = values[1], "soft_cos_sim" = values[2]
    ))
}

save(correlations_df, file = output_file_path)

