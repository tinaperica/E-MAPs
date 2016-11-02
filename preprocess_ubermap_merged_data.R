library(reshape2)

options(stringsAsFactors = F)

s.lim.point<-c(-3, 2)

setwd("~/Documents/GSP1/E_MAP_data/emap_analysis/")

### merged and averaged E-MAP data for Gsp1 mutants
e.map<-read.delim("avg_merged_June2016_screen_for_Gia.txt", head=T)   ## use export for Gia because it has ORF names for library
e.map.tab<-melt(e.map, id.vars=c("Gene"), variable.name="library", value.name="score")
e.map<-e.map.tab[order(e.map.tab$Gene, e.map.tab$library, decreasing=T),]
e.map.temp<-data.frame(lapply(e.map, gsub, pattern = "GSP1:", replacement = "", perl = T))
e.map<-cbind(e.map[,2:3], "Gene" = e.map.temp$Gene)
query <- unique(e.map$Gene)
### ORF to gene name annotation from SGD
orf_gene_name_index<-read.delim("orf_gene_GO_sgd_annotation.txt", head=F)
orf_index<-unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)
#############################
### Gsp1 mutants merged with ubermap (gene names) - this is only necessary to distinguish point mutants (they are all the same in the ORF file)
temp.ubermap<-read.delim("gene_names_merge_w_Ubermap_500.txt", head = T, skip = 0, nrow = length(query))
ubermap.gene_names<-melt(temp.ubermap, id.vars = c("Gene"), variable.name = "library.ORF", valuename = "score")
mutants<-as.character(unique(ubermap.gene_names$Gene))
mutants<-gsub(mutants, pattern = "GSP1 - ", replacement = "", perl = T)
rm(temp.ubermap)
### Gsp1 mutants merged with ubermap (orf names)
# replace Gsp1 ORFs with mutation names (because they are all the same orf!)
temp.ubermap<-read.delim("orf_names_merge_w_Ubermap_500.txt", head = T, stringsAsFactors = F)
#### "YDR113C" "YDL155W" "YGR108W" "YGR109C" "YPR119W" "YPR120C" -> these Gene names occur multiple times (2 to 4)
temp.ubermap[1:length(query),1] <- mutants
## mutants only
temp.ubermap.mutants.only <- temp.ubermap[1:length(query),]
## all but mutants
temp.ubermap.ubergenes.only <- temp.ubermap[(length(query)+1):length(temp.ubermap[,1]),]
rm(temp.ubermap)
#write.table(temp.ubermap.ubergenes.only, file = "ubergenes_to_cluster.txt", quote = F, row.names = F, sep = "\t")
#### sort out those few Genes that occur multiple times (add a unique identifier to them)
### first find those Genes that are not unique (probably data from separate experiments/labs for the same thing? Stll better not to mix them up)
control.count <- data.frame(table(temp.ubermap.ubergenes.only$Gene))
not_uniq <- as.character(control.count$Var1[control.count$Freq > 1])  # this is a vector of all genes that occur more than once
### now find their corresponding rows and add a unique prefix number to them
for (temp.gene in not_uniq) {
  row.names <- which(temp.ubermap.ubergenes.only$Gene == temp.gene)
  for (i in 1:length(row.names)) {
    temp.ubermap.ubergenes.only$Gene[row.names[i]] <- paste(i, temp.ubermap.ubergenes.only$Gene[row.names[i]], sep = "_")
  }
}
### now those few non-unique Genes look like this: 1_YDL155W - YDL155W
## melt and clean all but mutants
ubergenes.ubermap.orf_genes <- melt(temp.ubermap.ubergenes.only, id.vars = c("Gene"), variable.name = "library.ORF", value.name = "score")
#ubergenes.ubermap.orf_genes.complete <- ubergenes.ubermap.orf_genes[complete.cases(ubergenes.ubermap.orf_genes), ]
rm(temp.ubermap.ubergenes.only)
### this replacement is to remove the redundant ORF before the unique ORF (e.g. YFL008W_TSQ321 - YFL008W_TSQ321 -> keep only the first one) - > first one has the unique i_ for the few redundant ones!
ubergenes.ubermap.clean <- data.frame(lapply(ubergenes.ubermap.orf_genes[,1:2], gsub, pattern = " - .+$", replacement = "", perl = T))
### this is to replace the YOR298C.A with YOR298C-A (R puts . instead of - when loading library genes as column names)
ubergenes.ubermap.clean <- data.frame(lapply(ubergenes.ubermap.clean[,1:2], gsub, pattern = "\\.", replacement = "-", perl = T))
ubergenes.ubermap.orf_genes.uniq <- ubergenes.ubermap.clean  #### these Gene names and library names should be as short as possible, but still unique (e.g. if there are multiple TS mutants, of the same ORF, distinguish them)
#### in this step I only keep the basic ORF - this is necessary so I can math to SGD gene name
ubergenes.ubermap.clean <- data.frame(lapply(ubergenes.ubermap.clean, gsub, pattern = "^[0-9]+_", replacement = "", perl = T))  ### this removes the unique i_ for the redundant ones
ubergenes.ubermap.clean <- data.frame(lapply(ubergenes.ubermap.clean, gsub, pattern = "_.+$", replacement = "", perl = T))  ## it is important to keep this order of replacements of _!!!
ubergenes.ubermap <- cbind(ubergenes.ubermap.clean, "score" = ubergenes.ubermap.orf_genes$score,  "Gene_uniq" = ubergenes.ubermap.orf_genes.uniq$Gene)
names(ubergenes.ubermap)[1] <- "Gene.ORF"
ubergenes.ubermap.gene_names <- merge(ubergenes.ubermap, orf_index, by.x = "Gene.ORF", by.y = "orf")
names(ubergenes.ubermap.gene_names)[5] <- "Gene.gene_name"
ubergenes.ubermap.gene_names_both <- merge(ubergenes.ubermap.gene_names, orf_index, by.x = "library.ORF", by.y = "orf")
names(ubergenes.ubermap.gene_names_both)[6] <- "library.gene_name"
ubergenes.ubermap.gene_names_both <- ubergenes.ubermap.gene_names_both[order(ubergenes.ubermap.gene_names_both$Gene_uniq, ubergenes.ubermap.gene_names_both$library.ORF),]
write.table(ubergenes.ubermap.gene_names_both, file = "preprocessed_ubermap_ubergenes_only.txt", sep = "\t", quote = F, row.names = F)
# melt and clean mutants only
mut.ubermap <- melt(temp.ubermap.mutants.only, id.vars = c("Gene"), variable.name = "library.ORF", value.name = "score")
mut.ubermap.clean <- data.frame(lapply(mut.ubermap[,1:2], gsub, pattern = "\\.", replacement = "-", perl = T))
mut.ubermap <- cbind(mut.ubermap.clean, "score" = mut.ubermap$score)
rm(temp.ubermap.mutants.only)
mut.ubermap <- cbind(mut.ubermap, "ORF" = "YLR293C")
write.table(mut.ubermap, file = "preproceessed_ubermap_mut_only.txt", sep = "\t", quote = F, row.names = F)

#### make a subset of ubergenes.ubermap.gene_names_both that only has library genes that have scores outside the s.lim.point with at least one of the mutants
mut_ubermap_significant_scores <- mut.ubermap[ findInterval( mut.ubermap$score, s.lim.point ) != 1, ]
mut_significant_library_genes <- as.character(unique(mut_ubermap_significant_scores$library.ORF))  ## 1046 genes
all_library_genes <- as.character(unique(mut.ubermap$library.ORF)) #### 1356 genes

write.table(mut_ubermap_significant_scores, file = "preprocessed_ubermap_mut_only_significant.txt", sep = "\t", quote = F, row.names = F)

ubergenes_ubermap_significant_scores <- ubergenes.ubermap.gene_names_both[ ubergenes.ubermap.gene_names_both$library.ORF %in% mut_significant_library_genes, ]
write.table(ubergenes_ubermap_significant_scores, file = "preprocessed_ubermap_ubergenes_only_significant.txt", sep = "\t", quote = F, row.names = F)


discarded <- ubergenes.ubermap.gene_names_both[ ! ubergenes.ubermap.gene_names_both$library.ORF %in% mut_significant_library_genes, ]
discarded.library.genes <- as.character(unique(discarded$library.gene_name))

#### combine ubergenes and mutants only ubermap into all.ubermap before exporting
temp.mut.ubermap <- mut.ubermap
names(temp.mut.ubermap) <- c("Gene_uniq", "library", "score", "ORF")
temp.ubergenes.ubermap <- ubergenes.ubermap[c(4,2,3,1)]
names(temp.ubergenes.ubermap) <- c("Gene_uniq", "library", "score", "ORF")
all.ubermap <- rbind(temp.mut.ubermap, temp.ubergenes.ubermap)
##########################################################
write.table(all.ubermap, file = "preprocessed_ubermap_all.txt", sep = "\t", quote = F, row.names = F)


all_ubermap_significant_scores <- all.ubermap[ all.ubermap$library %in% mut_significant_library_genes, ]
write.table(all_ubermap_significant_scores, file = "preprocessed_ubermap_all_significant.txt", sep = "\t", quote = F, row.names = F)
