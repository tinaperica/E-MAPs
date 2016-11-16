### this script does all the clusters but for a small subset (56242) of pairs per task

library(getopt, lib.loc = "/netapp/home/tperica/R/x86_64-redhat-linux-gnu-library/3.3")

options(stringsAsFactors = F)

n.goods.cutoff<-5
s.lim.point<-c(-3, 2)
ratio.lim<-c(0.3, 3)

#### load an orf to gene name index
orf_gene_name_index <- read.delim("orf_gene_GO_sgd_annotation.txt", head = F)
orf_index <- unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)

#### library genes clustered
cluster_method <- "library_pearson_complete" #
cluster_file <- paste0("2016-10-21_", cluster_method, "_clusters.txt")
library_clusters <- read.delim(cluster_file, head = T)
clusters <- 1:max(library_clusters$cluster)
##############


## load preprocessed ubermap data (preprocessed so that ubermap$Gene is either a Gsp1 mutant or an ORF)
ubermap <- read.delim("preprocessed_ubermap_all.txt", head = T)

all_genes_and_mutants_df <- read.delim("genes_and_mutants_to_test.txt", head = F)
all_genes_and_mutants <- all_genes_and_mutants_df$V1

first_pair <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) )   #
last_pair <- as.numeric( Sys.getenv( "SGE_TASK_STEPSIZE" ) )

#library_clusters <- library_clusters[library_clusters$cluster %in% clusters,]
ubermap <- ubermap[ubermap$library %in% library_clusters$ORF,]
xy.limits <- range(ubermap$score, na.rm = T)

outputfilename <- paste0(Sys.Date(), "_", cluster_method, "_", first_pair, "_correlation_network_clusters.RData")
output_file_path <- file.path("Output", cluster_method, outputfilename)

pairs<-combn(all_genes_and_mutants, 2)

### for code testing purposes
#query1 <- "R108L"
#q1 = query1
#query2 <- "YJR074W"
#q2 <- query2
#query2 <- "YGL217C"
#query2 <- "YMR235C_TSQ172"
correlations_df <- data.frame()
for (i in seq_along(clusters)) {
  cluster = clusters[i]
  for ( p in first_pair:last_pair ) {
    query1 <- pairs[1, p]
    query2 <- pairs[2, p]
    temp_library_clusters <- library_clusters[library_clusters$cluster == cluster,]
    # grep, don't match, because I need to make sure to fit Gene_uniq
    ### e.g. cases like YMR235C_TSQ172
    ### also make sure to be able to handle cases when there are multiple mutants of the same gene - grep can have more than 1 hit
    temp.ubermap.query1 <- ubermap[ubermap$library %in% temp_library_clusters$ORF & grepl(query1, ubermap$Gene_uniq),]
    if (length(temp.ubermap.query1[,1]) > 0) {
      temp.ubermap.query2 <- ubermap[ubermap$library %in% temp_library_clusters$ORF & grepl(query2, ubermap$Gene_uniq),]
      if (length(temp.ubermap.query2) > 0) {
        #### queries1 and queries2 are vectors, usually of length 1, but can be more, if e.g. more than one TS mutant for gene
        queries1 <- as.character(unique(temp.ubermap.query1$Gene_uniq))
        queries2 <- as.character(unique(temp.ubermap.query2$Gene_uniq))
        for (q1 in queries1) {
          for (q2 in queries2) {
            temp.ubermap.q1 <- temp.ubermap.query1[temp.ubermap.query1$Gene_uniq == q1,]
            temp.ubermap.q2 <- temp.ubermap.query2[temp.ubermap.query2$Gene_uniq == q2,]
            merged.for.cor <- merge(temp.ubermap.q1, temp.ubermap.q2, by = "library")
            merged.for.cor <- merged.for.cor[complete.cases(merged.for.cor),]
            names(merged.for.cor) <- c("library", "Gene_uniq_q1", "score_q1", "ORF_q1", "Gene_uniq_q2", "score_q2", "ORF_q2")
            ### check if there is anything on the diagonal
            diagonal.merged.for.cor <- merged.for.cor[ findInterval(merged.for.cor$score_q1, s.lim.point) != 1 & 
                                                         findInterval(merged.for.cor$score_q2, s.lim.point) != 1, ]
            if ( length( diagonal.merged.for.cor[, 1]) > n.goods.cutoff ) {
              gene_names_q1 <- merge(orf_index, merged.for.cor, by.x = "orf", by.y = "ORF_q1" )
              gene_names_q1 <- within(gene_names_q1, "q1" <- paste(gene_name, Gene_uniq_q1, sep = " - "))
              if (length(gene_names_q1[,1]) > 0) {
                gene_names_q1_q2 <- merge(orf_index, gene_names_q1[, c(3,5,6,7,8,9)], by.x = "orf", by.y = "ORF_q2")
                gene_names_q1_q2 <- within(gene_names_q1_q2, "q2" <- paste(gene_name, Gene_uniq_q2, sep = " - "))
                gene_names_q1_q2 <- gene_names_q1_q2[, c(3, 4, 6, 7, 8)]
                if (length(gene_names_q1_q2[,1]) > 0) {
                  correlation <- cor(gene_names_q1_q2$score_q1, gene_names_q1_q2$score_q2, use="pairwise.complete.obs")
                  gene_names_q1_q2 <- gene_names_q1_q2[ findInterval(gene_names_q1_q2$score_q1, s.lim.point) != 1 & 
                                                          findInterval(gene_names_q1_q2$score_q2, s.lim.point) != 1, ]
                  gene_names_q1_q2 <- cbind(gene_names_q1_q2, "score_ratio" = gene_names_q1_q2$score_q1/gene_names_q1_q2$score_q2)
                  correlations_df <- rbind(
                    correlations_df, 
                    data.frame( 
                      "gene1" = gene_names_q1_q2$q1[1], "gene_2" = gene_names_q1_q2$q2[1], 
                      "cluster_method" = cluster_method, "cluster" = cluster,
                      "corr"  = round(correlation, 3)
                    )
                  )
                }
              }
            }
          } 
        }
      }
    }
  }
}
write.table(correlations_df, file = output_file_path, quote = F, sep = "\t", row.names = F)
