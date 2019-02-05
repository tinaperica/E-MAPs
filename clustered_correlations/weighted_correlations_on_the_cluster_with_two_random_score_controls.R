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

load(inputfilename)
load("20180630_emap_data_for_corr.RData")
outputfilename <- paste0("20180702_", task_n, "_correlations.RData")
output_file_path <- file.path("Output/correlations", outputfilename)
correlations_df <- data.frame()
for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  ubermap.query1 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query1, ]
  ubermap.query2 <- ubermap[["ubermap"]][ubermap[["ubermap"]][["Gene_uniq"]] == query2, ]
  for ( i in seq_along( ubermap[["clusters"]] ) ) {
    cluster <- ubermap[["clusters"]][i]
    temp_library_clusters <- ubermap[["library_clusters"]][ubermap[["library_clusters"]][["cluster"]] == cluster, ]
    cluster_size <- length(temp_library_clusters$ORF)
    temp.ubermap.query1 <- ubermap.query1[ubermap.query1[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    temp.ubermap.query2 <- ubermap.query2[ubermap.query2[["library.ORF"]] %in% temp_library_clusters[["ORF"]], ]
    merged.emap.data <- merge(temp.ubermap.query1[, c("Gene_uniq", "library.ORF","score","random_score", "random_score_2", "random_high_score", "random_high_score_2")], 
                              temp.ubermap.query2[, c("Gene_uniq", "library.ORF", "score", "random_score", "random_score_2", "random_high_score", "random_high_score_2")], 
                              by = "library.ORF")
    data_completness <- merged.emap.data[, c("score.x", "score.y")]
    data_completness <- data_completness[complete.cases(data_completness),]
    data_completness <- length(data_completness$score.x) / cluster_size
    if (data_completness > 0.5) {
      values <- round(compare_wrap(merged.emap.data$score.x, merged.emap.data$score.y), 4)
      random_values <- round(compare_wrap(merged.emap.data$random_score.x, merged.emap.data$random_score.y), 4)
      random_values_2 <- round(compare_wrap(merged.emap.data$random_score_2.x, merged.emap.data$random_score_2.y), 4)
      random_high_values <- round(compare_wrap(merged.emap.data$random_high_score.x, merged.emap.data$random_high_score.y), 4)
      random_high_values_2 <- round(compare_wrap(merged.emap.data$random_high_score_2.x, merged.emap.data$random_high_score_2.y), 4)
      correlations_df <- rbind(correlations_df, data.frame( 
          "Gene_uniq1" = query1, "Gene_uniq2" = query2, "cluster" = cluster, 
          "w_correlation" = values[1], "soft_cos_sim" = values[2], 
          "random_w_correlation" = mean(c(random_values[1], random_values_2[1]), na.rm = T),
          "random_soft_cos_sim" = mean(c(random_values[2], random_values_2[2]), na.rm = T),
          "random_high_w_correlation" = mean(c(random_high_values[1], random_high_values_2[1]), na.rm = T),
          "random_high_soft_cos_sim" = mean(c(random_high_values[2], random_high_values_2[2]), na.rm = T)
        ))
    }
  }
}

save(correlations_df, file = output_file_path)

