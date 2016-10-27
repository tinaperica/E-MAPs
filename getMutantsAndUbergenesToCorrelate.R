### Ubermap correlations 
#### run on guybrush in the E-MAP/June2016_analysis dir
#### screen -S RCorPairs
#### Rscript getMutantsAndUbergenesToCorrelate.R
library(reshape2)
library(calibrate)
library(plyr)
setwd("/home/tina/E-MAP/June2016_analysis")
#setwd("~/Documents/GSP1/E_MAP_data/June2016_analysis/")
orf_gene_name_index<-read.delim("orf_gene_GO_sgd_annotation.txt", head=F)
orf_index<-unique(data.frame("orf" = orf_gene_name_index$V1, "gene_name" = orf_gene_name_index$V2))
rm(orf_gene_name_index)
e.map<-read.delim("avg_merged_June2016_screen_for_Gia.txt", head=T, sep="\t")   ## use export for Gia because it has ORF names for library
e.map.tab<-melt(e.map, id.vars=c("Gene"), variable.name="library", value.name="score")
e.map<-e.map.tab[order(e.map.tab$Gene, e.map.tab$library, decreasing=T),]
e.map.temp<-data.frame(lapply(e.map, gsub, pattern = "GSP1:", replacement = "", perl = T))
e.map<-cbind(e.map[,2:3], "Gene" = e.map.temp$Gene)
query<-levels(e.map$Gene)
lib.genes<-levels(e.map$library)
#slope.cutoff<-0.1
n.goods.cutoff<-5
s.lim.point<-c(-3, 2)  ### this is also the threshold for the wt control
s.lim.text<-c(-5, 4)
ratio.lim<-c(0.3, 3)
x.limits=c(min(e.map$score, na.rm=TRUE), max(e.map$score, na.rm=TRUE))
y.limits=c(min(e.map$score, na.rm=TRUE), max(e.map$score, na.rm=TRUE))
ubermap.merged.gene_name_file<-"gene_names_merge_w_Ubermap_500.txt"
temp.ubermap<-read.delim(ubermap.merged.gene_name_file, head = T, skip = 0, nrow = length(query), stringsAsFactors = F)
ubermap.gene_names<-melt(temp.ubermap, id.vars = c("Gene"), variable.name = "library", valuename = "score")
rm(temp.ubermap)
mutants<-as.character(unique(ubermap.gene_names$Gene))
mutants<-gsub(mutants, pattern = "GSP1 - ", replacement = "", perl = T)
ubermap.merged.orf_file<-"orf_names_merge_w_Ubermap_500.txt"
temp.ubermap<-read.delim(ubermap.merged.orf_file, head = T, stringsAsFactors = F)
temp.ubermap[1:length(query),1] <- mutants
ubermap.orf_names<-melt(temp.ubermap, id.vars=c("Gene"), variable.name="library", value.name="score")
ubermap.orf_names.complete <- ubermap.orf_names[complete.cases(ubermap.orf_names), ]
rm(temp.ubermap)
ubermap.clean<-data.frame(lapply(ubermap.orf_names.complete[,1:2], gsub, pattern = " - .+$", replacement = "", perl = T))
ubermap.clean<-data.frame(lapply(ubermap.clean, gsub, pattern = "_.+$", replacement = "", perl = T))
ubermap.clean<-cbind(ubermap.clean, "score" = ubermap.orf_names.complete$score)
e.map.wt<-subset(e.map, (Gene == "GSP1-NAT" & findInterval(score, s.lim.point) == 1))
genes_and_mutants<-unique(as.character(ubermap.clean$Gene))
genes <- genes_and_mutants[! genes_and_mutants %in% mutants]
genes_and_mutants_to_test <- mutants
for (m in 1:length(mutants)) {
  mut <- mutants[m]
  filename = paste(mut, "_ubermap_pairwise_corr.pdf", collapse="")
  pdf(file = filename)
  ubermap.mut <- subset(ubermap.clean, Gene == mut)
  ubermap.mut <- merge(ubermap.mut, e.map.wt, by = "library")
  ubermap.mut <- ubermap.mut[,1:3]
  names(ubermap.mut) <- c("library", "Gene", "score")
  for (ug in 1:length(genes)) {
    gene <- genes[ug]
    #gene = "YGR119C"
    ubermap.ubergene <- subset(ubermap.clean, Gene == gene)
    ubermap.ubergene <- merge(ubermap.ubergene, e.map.wt, by = "library")
    ubermap.ubergene <- ubermap.ubergene[,1:3]
    names(ubermap.ubergene) <- c("library", "Gene", "score")
    ubermap.final.merge <- merge(ubermap.mut, ubermap.ubergene, by = "library")
    if (length(ubermap.final.merge$score.x) > 500) {
      correlations.table <- ddply(ubermap.final.merge, "Gene.x", function(df) cor(df$score.x, df$score.y, use = "pairwise.complete.obs"))
      pvalue.table <- ddply(ubermap.final.merge, "Gene.x", function(df) cor.test(df$score.x, df$score.y)$p.value)
      labels.x <- ubermap.final.merge$score.x[abs(ubermap.final.merge$score.x) > 5 | abs(ubermap.final.merge$score.y) > 5]
      labels.y <- ubermap.final.merge$score.y[abs(ubermap.final.merge$score.x) > 5 | abs(ubermap.final.merge$score.y) > 5]
      labels.labels <- ubermap.final.merge$library[abs(ubermap.final.merge$score.x) > 5 | abs(ubermap.final.merge$score.y) > 5]
      plot(ubermap.final.merge$score.x, ubermap.final.merge$score.y, xlab = mut, ylab = gene)
      if (length(labels.labels) > 0) {text(labels.x, labels.y, labels=labels.labels, cex=0.6, pos=4)}
      if (correlations.table$V1 > 0.1 & pvalue.table$V1 < 0.1) {
        score.ratio <- with (ubermap.final.merge, abs(score.y/score.x))
        final <- cbind(ubermap.final.merge, score.ratio)
        diagonal <- subset(final, findInterval(score.ratio, ratio.lim) == 1 & findInterval(score.y, s.lim.point) != 1 & findInterval(score.x, s.lim.point) != 1)
        if (length(diagonal$library) > n.goods.cutoff) {
          genes_and_mutants_to_test <- append(genes_and_mutants_to_test, gene)
        }   
      }
    }
  }
  dev.off()
}
genes_and_mutants_to_test <- data.frame(unique(genes_and_mutants_to_test))
write.table(genes_and_mutants_to_test, file = "genes_and_mutants_to_test.txt", quote = F, row.names = F, sep = "\t")

