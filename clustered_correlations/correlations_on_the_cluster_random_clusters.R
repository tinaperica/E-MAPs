###### MODIFIED on 20180108 - It never calculates the correlation if fewer than 5 points (2 points will always have correlation of 1 or -1)
##### modified on 20180123 - the part where I remove the library genes that have significant correlations is now moved to the preorocessing script
############# no point to keep doing that for each pair of genes separately
### this script does all the clusters but for a small subset of pairs per task (times 35 clusters)

options(stringsAsFactors = F)

load("20180122_random_control_correlations_pair_preload.RData")

first_pair <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
step_size <- as.numeric( Sys.getenv( "SGE_TASK_STEPSIZE" ) ) - 1
last_pair <- first_pair + step_size
outputfilename <- paste0("20180131_", ubermap[["cluster_method"]], "_", first_pair, "_correlation_network_clusters.RData")
output_file_path <- file.path("Output", ubermap[["cluster_method"]], outputfilename)

orfs_to_plot <- c(
  "YJR074W", "YDR002W", "YMR235C", "YGL097W", "YMR308C"
)

correlations_df <- data.frame()

for ( p in first_pair:last_pair ) {
  query1 <- ubermap[["pairs"]][1, p]
  query2 <- ubermap[["pairs"]][2, p]
  
  # grep, don't match, because I need to make sure to fit Gene_uniq
  ### e.g. cases like YMR235C_TSQ172
  ### also make sure to be able to handle cases when there are multiple mutants of the same gene - grep can have more than 1 hit
  ubermap.query1 <- ubermap[["ubermap"]][grepl( query1, ubermap[["ubermap"]][["Gene_uniq"]] ), ]
  ubermap.query2 <- ubermap[["ubermap"]][grepl( query2, ubermap[["ubermap"]][["Gene_uniq"]] ), ]
  
  if ( (length(orfs_to_plot[grepl(query1, orfs_to_plot)]) > 0 | length(orfs_to_plot[grepl(query2, orfs_to_plot)]) > 0) &
       (ubermap.query1$ORF[1] == "YLR293C" | ubermap.query1$ORF[1] == "YLR293C") ) {
    plot <- "yes"
    pdf(file = paste0("pdf_output/", ubermap[["cluster_method"]], "_", query1, "_", query2, "_scatterplots.pdf"))
  } else {
    plot <- "no"
  }
  
  for ( i in seq_along( ubermap[["clusters"]] ) ) {
    
    cluster <- ubermap[["clusters"]][i]
    temp_library_clusters <- ubermap[["library_clusters"]][ubermap[["library_clusters"]][["cluster"]] == cluster, ]
    temp.ubermap.query1 <- ubermap.query1[ubermap.query1[["library"]] %in% temp_library_clusters[["ORF"]], ]
    
    if (length(temp.ubermap.query1[,1]) > 5) {   #### I changed all of these conditions to > 5, as that is the final condition below, no point in keeping this as just > 0
      
      temp.ubermap.query2 <- ubermap.query2[ubermap.query2[["library"]] %in% temp_library_clusters[["ORF"]], ]
      
      if (length(temp.ubermap.query2[,1]) > 5) {
        
        #### queries1 and queries2 are vectors, usually of length 1, but can be more, if e.g. more than one TS mutant for gene
        queries1 <- as.character(unique(temp.ubermap.query1$Gene_uniq))
        queries2 <- as.character(unique(temp.ubermap.query2$Gene_uniq))
        
        for (q1 in queries1) {
          for (q2 in queries2) {
            
            temp.ubermap.q1 <- temp.ubermap.query1[temp.ubermap.query1$Gene_uniq == q1,]
            temp.ubermap.q2 <- temp.ubermap.query2[temp.ubermap.query2$Gene_uniq == q2,]
            gene_names_q1 <- merge(ubermap[["orf_index"]], temp.ubermap.q1, by.x = "orf", by.y = "ORF" )
            gene_names_q1 <- within(gene_names_q1, "q1" <- paste(gene_name, Gene_uniq, sep = " - "))
            
            if (length(gene_names_q1[,1]) > 5) {
              
              gene_names_q2 <- merge(ubermap[["orf_index"]], temp.ubermap.q2, by.x = "orf", by.y = "ORF" )
              gene_names_q2 <- within(gene_names_q2, "q2" <- paste(gene_name, Gene_uniq, sep = " - "))
              gene_names_q1_q2 <- merge(gene_names_q1[, c(4,5,6)], gene_names_q2[, c(4,5,6)], by = "library")
              gene_names_q1_q2_complete_cases <- gene_names_q1_q2[complete.cases(gene_names_q1_q2), ]
              
              if (length(gene_names_q1_q2_complete_cases[,1]) > 5) { # on 20180124 this was changed from 0 to 5 - because if you have only 2 points, corr will always be 1!
                
                names(gene_names_q1_q2_complete_cases) <- c("library", "score_q1", "q1", "score_q2", "q2")
                if (plot == "yes") {
                  plot(gene_names_q1_q2_complete_cases$score_q1, gene_names_q1_q2_complete_cases$score_q2, main = cluster,
                       xlab = gene_names_q1_q2_complete_cases$q1[1], ylab = gene_names_q1_q2_complete_cases$q2[1])
                }
                
                gene_names_q1_q2_complete_cases <- cbind( gene_names_q1_q2_complete_cases, "score_ratio" = abs(gene_names_q1_q2_complete_cases$score_q1)/abs(gene_names_q1_q2_complete_cases$score_q2) )
                correlation <- cor(gene_names_q1_q2_complete_cases$score_q1, gene_names_q1_q2_complete_cases$score_q2, use="pairwise.complete.obs")
                total.n.genes <- length(gene_names_q1_q2_complete_cases[complete.cases(gene_names_q1_q2_complete_cases),][,1])
                correlations_df <- rbind(
                  correlations_df, 
                  data.frame( 
                    "gene1" = gene_names_q1_q2_complete_cases$q1[1], "gene2" = gene_names_q1_q2_complete_cases$q2[1], 
                    "cluster_method" = ubermap[["cluster_method"]], "cluster" = cluster, "corr"  = round(correlation, 3) )
                )
              } else {
                correlations_df <- rbind(
                  correlations_df,
                  data.frame(
                    "gene1" = gene_names_q1_q2$q1[1], "gene2" = gene_names_q1_q2$q2[1],
                    "cluster_method" = ubermap[["cluster_method"]], "cluster" = cluster, "corr"  = "NA")
                )
              }
            }
          } 
        }
      }
    }
  }
  if (plot == "yes") {
    dev.off()
  }
}
save(correlations_df, file = output_file_path)
