library(gplots)

options(stringsAsFactors = F)

setwd("~/Documents/GSP1/E_MAP_data/emap_analysis/clustered_correlations/correlations_data/RData")

ms_partners_df <- read.delim("../../../Gsp1_MS_hits.txt", head = F)
ms_partners <- ms_partners_df$V1


dir_list <- dir()

all_correlations <- data.frame()

for (d in seq_along(dir_list)) {
  dir <- dir_list[d]
  files_list <- list.files(path = dir)
  
  for (f in seq_along(files_list)) {
    filename <- files_list[f]
    if (grepl(".RData", filename)) {
      load(file.path(dir, filename))
      all_correlations <- rbind(all_correlations, correlations_df)
    }
  }
  if (length(all_correlations[1,]) > 0) {
    rm(correlations_df)
    names(all_correlations) <- c("gene1", "gene2", "clustering_method", "cluster", "corr")
    write.table(all_correlations, file = paste0(dir, "/", Sys.Date(), "_merged_correlations.txt"), quote = F, sep = "\t", row.names = F)
  }
}
  all_genes <- as.character(unique(c(all_correlations$gene1, all_correlations$gene2)))
  mutants <- all_genes[grepl("GSP1", all_genes)]
  partners <- all_genes[! all_genes %in% mutants]
  mutants <- gsub(mutants, pattern = "GSP1 - ", replacement = "", perl = T)
  partners <- gsub(partners, pattern = " - .+$", replacement = "", perl = T)
  
  clean_names <- data.frame(lapply(all_correlations[,1:2], gsub, pattern = "GSP1 - ", replacement = "", perl = T))
  clean_names <- data.frame(lapply(clean_names, gsub, pattern = " - .+$", replacement = "", perl = T))
  
  correlations <- cbind(clean_names, all_correlations[,3:5])
  clustering_method <- dir
  
  
  mutant_partner_corr_of_corr <- data.frame()
  clusters <- as.character(unique(correlations$cluster))
  for (m in seq_along(mutants)) {
    mutant = mutants[m]
    GI_mut <- data.frame()
    GI_mut1 <- correlations[ correlations$gene1 == mutant, ]
    GI_mut1 <- GI_mut1[ ! GI_mut1$gene2 %in% mutants, ]
    GI_mut <- rbind(GI_mut, data.frame("gene" = GI_mut1$gene2, "cluster" = GI_mut1$cluster,  "gene_cluster" = paste(GI_mut1$gene2, GI_mut1$cluster), "corr" = GI_mut1$corr))
    GI_mut2 <- correlations[ correlations$gene2 == mutant, ]
    GI_mut <- rbind(GI_mut, data.frame("gene" = GI_mut2$gene1, "cluster" = GI_mut2$cluster, "gene_cluster" = paste(GI_mut2$gene1, GI_mut2$cluster), "corr" = GI_mut2$corr))
    GI_mut <- GI_mut[order(as.character(GI_mut$gene, GI_mut$cluster)),]
    ### remove correlations with other mutants
    GI_mut <- GI_mut[ ! GI_mut$gene %in% mutants, ]
    for ( p in seq_along(partners) ) {
      partner = partners[p]
      ### remove the correlation with that partner from the GI_mut
      GI_mut_temp <- GI_mut[ ! GI_mut$gene %in% partner, ]
      GI_partner <- data.frame()
      GI_partner1 <- correlations[ correlations$gene1 == partner, ]
      GI_partner2 <- correlations[ correlations$gene2 == partner, ]
      if (length(GI_partner1[,1]) > 0 | length(GI_partner2[,1]) > 0) {
        GI_partner <- rbind(GI_partner, data.frame("gene" = GI_partner1$gene2, "cluster" = GI_partner1$cluster, 
                                                   "gene_cluster" = paste(GI_partner1$gene2, GI_partner1$cluster), "corr" = GI_partner1$corr))
        GI_partner <- rbind(GI_partner, data.frame("gene" = GI_partner2$gene1, "cluster" = GI_partner2$cluster,
                                                   "gene_cluster" = paste(GI_partner2$gene1, GI_partner2$cluster), "corr" = GI_partner2$corr))
        ### remove correlations with other mutants
        GI_partner <- GI_partner[ ! GI_partner$gene %in% mutants, ]
        partner_specific <- setdiff(GI_partner$gene_cluster, GI_mut_temp$gene_cluster)
        GI_partner <- GI_partner[ ! GI_partner$gene_cluster %in% partner_specific, ]
        mutant_specific <- setdiff(GI_mut_temp$gene_cluster, GI_partner$gene_cluster)
        GI_mut_temp <- GI_mut_temp[ ! GI_mut_temp$gene_cluster %in% mutant_specific, ]
        GI_partner <- GI_partner[ order( as.character(GI_partner$gene_cluster) ), ]
        GI_mut_temp <- GI_mut_temp[ order( as.character(GI_mut_temp$gene_cluster) ), ]
        temp.measure.df <- data.frame()
        partner.sd <- sd(GI_partner$corr, na.rm = T)
        mutant.sd <- sd(GI_mut_temp$corr, na.rm = T)
        corr_of_corr <- 0
        if (abs(partner.sd) > 0 & abs(mutant.sd) > 0) {
          temp_for_corr <- merge(GI_partner, GI_mut_temp, by = "gene_cluster")
          corr_of_corr <- cor(temp_for_corr$corr.x, temp_for_corr$corr.y, method = "pearson", use = "pairwise.complete.obs")
        }
        temp.measure.df <- rbind(temp.measure.df, data.frame(mutant, partner, "measure" = "corr_of_corr", "cluster" = "all", "value" = corr_of_corr))
        
        for (clust in clustßers) {
          GI_partner_clust <- GI_partner[GI_partner$cluster == clust, ]
          GI_mut_clust <- GI_mut_temp[GI_mut_temp$cluster == clust, ]
          
          partner_specific <ß- setdiff(GI_partner_clust$gene_cluster, GI_mut_clust$gene_cluster)
          GI_partner_clust <- GI_partner_clust[ ! GI_partner_clust$gene_cluster %in% partner_specific, ]
          mutant_specific <- setdiff(GI_mut_clust$gene_cluster, GI_partner_clust$gene_cluster)
          GI_mut_clust <- GI_mut_clust[ ! GI_mut_clust$gene_cluster %in% mutant_specific, ]
          GI_partner_clust <- GI_partner_clust[ order( as.character(GI_partner_clust$gene_cluster) ), ]
          GI_mut_clust <- GI_mut_clust[ order( as.character(GI_mut_clust$gene_cluster) ), ]
          
          partner.sd <- sd(GI_partner_clust$corr, na.rm = T)
          mutant.sd <- sd(GI_mut_clust$corr, na.rm = T)
          corr_of_corr <- 0
          if (abs(partner.sd) > 0 & abs(mutant.sd) > 0) {
            temp_for_corr <- merge(GI_partner_clust, GI_mut_clust, by = "gene_cluster")
            corr_of_corr <- cor(temp_for_corr$corr.x, temp_for_corr$corr.y, method = "pearson", use = "pairwise.complete.obs")
          }
          temp.measure.df <- rbind(temp.measure.df, data.frame(mutant, partner, "measure" = "corr_of_corr", "cluster" = clust, "value" = corr_of_corr))
        }
        
        mutant_partner_corr_of_corr <- rbind(mutant_partner_corr_of_corr, temp.measure.df)
      }
    }
  }
  
  write.table(mutant_partner_corr_of_corr, file = paste0(clustering_method, "_corr_of_corr.txt"), quote = F, sep = "\t", row.names = F)
  
  clusters <- as.character(unique(mutant_partner_corr_of_corr$cluster))
  partners.list<- list()
  partners.list[["Partners_with_structures"]] <- c("RNA1", "SRM1", "YRB1", "PSE1", "MOG1")
  partners.list[["MS_partners"]] <- ms_partners
  #partners.list[["All_Partners"]] <- partners
  selected_measures <- c("corr_of_corr")
  corr.mutants <- as.character(unique(mutant_partner_corr_of_corr$mutant))
  corr.partners <- as.character(unique(mutant_partner_corr_of_corr$partner))
  for ( l in seq_along( partners.list ) ) {
    for ( m in selected_measures ) {
      for ( clust in seq_along(clusters) ) {
        cluster <- clusters[clust]
        temp_measure_df <- mutant_partner_corr_of_corr[mutant_partner_corr_of_corr$measure == m & mutant_partner_corr_of_corr$cluster == cluster, ]
        temp.df <- temp_measure_df[temp_measure_df$partner %in% unlist(partners.list[[l]]),]
        temp.df <- temp.df[ order(temp.df$mutant, temp.df$partner), ]
        temp.partners <- as.character( unique( temp.df$partner ))
        colnames <- as.character( unique( temp.df$partner ) ) 
        correlation.matrix <- matrix( temp.df$value, byrow = T, nrow = length(corr.mutants), ncol = length(temp.partners) )
        correlation.matrix <- t(correlation.matrix)
        colnames(correlation.matrix)<-corr.mutants
        rownames(correlation.matrix)<-colnames
        filename <- paste0(clustering_method, "_", names(partners.list[l]), "_", m, "_", cluster, ".pdf")
        pdf(file = filename)
        heatmap <- heatmap.2( correlation.matrix, scale = "none", dendrogram="both", 
                              trace="none", density.info="none", col = cm.colors, 
                              key.ylab="partners", key.xlab = m, 
                              main = names(partners.list[l]), cexCol=0.6, cexRow = 0.6, keysize=0.9, margin=c(6,6))
        dev.off()
      }
    }
  }
  
}

