library(gplots)

options(stringsAsFactors = F)

ms_partners_df <- read.delim("Gsp1_MS_hits.txt", head = F)
ms_partners <- ms_partners_df$V1

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
}
