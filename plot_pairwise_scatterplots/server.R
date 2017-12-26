#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)
colors <- c(brewer.pal(8, "Set1"), brewer.pal(11, "Spectral"))
ubermap <- readRDS("data/ubermap.rds")
library_clusters <- readRDS("data/library_clusters.rds")
clusters <- as.character(unique(library_clusters[["cluster"]]))
cluster_GO_df <- readRDS("data/cluster_GO_df.rds")
GO_clusters <- as.character(unique(cluster_GO_df[["GO_cluster"]]))
GO_cluster_cols <- list(GO_clusters, colors[1:length(GO_clusters)])
orf_gene_name_index <- read.delim("data/orf_gene_GO_sgd_annotation.txt", head = F)
names(orf_gene_name_index) <- c("orf", "gene_name", "GO")
uniq_ubermap_gene_names <- readRDS("data/unique_gene_names_from_ubermap.rds")

shinyServer( function(input, output) {
  
  output$scatterplot <- renderPlot( {
    input1 <- toupper(input$gene1)
    if (length(orf_gene_name_index$orf[orf_gene_name_index$gene_name == input1]) > 0 ) {
      orf1 <- orf_gene_name_index$orf[orf_gene_name_index$gene_name == input1][1] 
      gene1 <- input1
    } else if (length(orf_gene_name_index$gene_name[orf_gene_name_index$orf == input1]) > 0) {
      orf1 <- input1
      gene1 <- orf_gene_name_index$gene_name[orf_gene_name_index$orf == input1][1]
    } else {  ### this is the case for mutant
      orf1 <- input1
      gene1 <- input1
    }
    
    input2 <- toupper(input$gene2)
    if (length(orf_gene_name_index$orf[orf_gene_name_index$gene_name == input2]) > 0 ) {
      orf2 <- orf_gene_name_index$orf[orf_gene_name_index$gene_name == input2][1] 
      gene2 <- input2
    } else if (length(orf_gene_name_index$gene_name[orf_gene_name_index$orf == input2]) > 0) {
      orf2 <- input2
      gene2 <- orf_gene_name_index$gene_name[orf_gene_name_index$orf == input2][1]
    } else {  ### this is the case for mutant
      orf2 <- input2
      gene2 <- input2
    }
    
    uniq_gene1 <- uniq_ubermap_gene_names[grepl(orf1, uniq_ubermap_gene_names)]
    uniq_gene2 <- uniq_ubermap_gene_names[grepl(orf2, uniq_ubermap_gene_names)]
    #### get mfrow from number of plots 
    op_vector <- n2mfrow(length(uniq_gene1) * length(uniq_gene2) + 1) ### +1 is for the legend plot
    ### this n2mfrow puts more rows than columns, so switch it around
    op <- par(mfrow = c(op_vector[2], op_vector[1]))
    
    for (i in seq_along(uniq_gene1)) {
      for (j in seq_along(uniq_gene2)) {
        uniq_orf1 <- uniq_gene1[i]
        uniq_orf2 <- uniq_gene2[j]
        gene1_ubermap <- ubermap[ubermap[["Gene_uniq"]] == uniq_orf1, ][,1:3]
        gene2_ubermap <- ubermap[ubermap[["Gene_uniq"]] == uniq_orf2, ][,1:3]
        merge_genes <- merge(gene1_ubermap, gene2_ubermap, by = "library")
        names(merge_genes) <- c("library", "gene1", "score1", "gene2", "score2")
        merge_to_plot <- merge(merge_genes, library_clusters, by.x = "library", by.y = "ORF") 
        merge_to_plot <- merge(merge_to_plot, cluster_GO_df, by = "cluster")
        #merge_to_plot <- cbind(merge_to_plot, "col" = "red")
        plot(merge_to_plot[["score1"]], merge_to_plot[["score2"]], 
             xlab = paste0(gene1, " (", uniq_orf1,")"), ylab = paste0(gene2, " (", uniq_orf2,")"), 
             xlim = c(-20,6), ylim = c(-20, 6), pch = 19, col = "gray")
        GO_clusters_to_plot <- input$checkGO
        for (g in seq_along(GO_clusters)) {
          GO_cluster <- GO_clusters[g]
          if (GO_cluster %in% GO_clusters_to_plot) {
            temp <- merge_to_plot[merge_to_plot[["GO_cluster"]] == GO_cluster, ]
            points(temp$score1, temp$score2, col = GO_cluster_cols[[2]][g], pch = 19)
          }
        }
        plot(merge_to_plot[["score1"]], merge_to_plot[["score2"]], 
             xlab = "", ylab = "", xlim = c(-15,5), ylim = c(-15, 5), type = "n", axes =F)
        legend("center", legend = GO_cluster_cols[[1]], fill = GO_cluster_cols[[2]], bty = "n")
        
      }
    }
  
  })
 
})
