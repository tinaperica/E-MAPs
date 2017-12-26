library(reshape2)
library(RColorBrewer)
mean.na<-function (x) {
  mean(x, na.rm=T)
}
min.na<-function (x) {
  min(x, na.rm=T)
}
abs.max <- function(x) {
  max(abs(x), na.rm=T)
}
abs.min <- function(x) {
  min(abs(x), na.rm=T)
}
### complex annotation  - cluster the point mutations
load("basic_E-MAP_data/June2016_Gsp1_E-MAP_data.RData")
#gia.results<-read.delim("RESULTS_GTPases/avg_merged_June2016_screen_for_Gia_gia_all.txt", head=T)
gia.results<-read.delim("RESULTS_scores_uncurated_complex_annotation/avg_merged_June2016_screen_for_Gia_gia_all.txt", head=T)
#complex.annotation<-read.delim("tina_bioprocess_annotation.txt", head=T)
#gia.results<-read.delim("RESULTS_scores_tina_bioprocess_annotation/avg_merged_June2016_screen_for_Gia_gia_all.txt", head=T)
#complex.annotation<-read.delim("GO_slim_annotation.txt", head=T)
#gia.results<-read.delim("RESULTS_scores_GO_slims_annotation/avg_merged_June2016_screen_for_Gia_gia_all.txt", head=T)
#head(gia.results)
#filename = "GTPases_boxplots_with_msd_2d.pdf"
filename="all_complexes_boxplots_with_msd_direct_projections.pdf"
#filename="tina_bioprocesses_boxplots_with_msd.pdf"
#filename="GO_slims_boxplots_with_msd.pdf"
gia.p.value<-with(gia.results, aggregate(P.value, by=list(complex=Group.B), min))
gia.p.value.sort<-gia.p.value[order(gia.p.value$x, decreasing=F),]
complexes<-as.vector(gia.p.value.sort$complex)
e.map.tab<-melt(e.map, id.vars=c("Gene"), variable.name="library", value.name="score")
e.map.tab<-subset(e.map.tab, (Gene!="REVERTANT1" & Gene!="REVERTANT2" & Gene!="REVERTANT3"))
e.map<-e.map.tab[order(e.map.tab$Gene, e.map.tab$library, decreasing=T),]
wt.emap<-subset(e.map, Gene == "GSP1-NAT")
emap.min<-floor(min(e.map$score, na.rm=T))
emap.max<-ceiling(max(e.map$score, na.rm=T))
genes<-as.vector(unique((e.map$Gene)))
genes<-sort(genes, decreasing=F)
library.genes<-levels(e.map$library)
e.map.complex<-merge(e.map, complex.annotation, by="library")
wt.emap.complex<-merge(wt.emap, complex.annotation, by="library")
colors = c(brewer.pal(8, "Dark2"),brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"),brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"),brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired"))
library.genes.in.complexes<-vector()
pdf(file=filename, width=10)
for (i in 1:length(complexes)) {
  compl<-complexes[i]
  temp.emap<-subset(e.map.complex, complex==compl)
  temp.wt.emap<-subset(wt.emap.complex, complex==compl)
  temp.emap<-temp.emap[order(temp.emap$gene_name, temp.emap$Gene),]
  temp.emap.matrix<-matrix(temp.emap$score, nrow=length(genes))
  rownames(temp.emap.matrix)<-genes
  temp.emap.matrix[is.na(temp.emap.matrix)]<-0
  temp.emap.matrix<-scale(temp.emap.matrix)
  distance<-dist(temp.emap.matrix, method="euclidean")
  fit<-cmdscale(distance, eig=T, k=2)
  multiscale.gene.order.table<-data.frame(genes, data.frame(dim1.order.scale=as.vector(fit$points[,1]), dim2.order.scale=as.vector(fit$points[,2])))
  multiscale.gene.dim1.order.table<-multiscale.gene.order.table[order(multiscale.gene.order.table$dim1.order.scale),]
  multiscale.gene.dim2.order.table<-multiscale.gene.order.table[order(multiscale.gene.order.table$dim2.order.scale, decreasing=T),]
  ordered_genes.dim1<-as.vector(multiscale.gene.dim1.order.table$genes)
  ordered_genes.dim2<-as.vector(multiscale.gene.dim2.order.table$genes)
  by.score<-with(temp.emap, aggregate(score, by=list(gene_name), abs.max))
  low.score.library.genes<-vector()
  if (length(by.score$Group.1[by.score$x < 4]) > 0) {
    low.score.library.genes<-as.vector(unique((by.score$Group.1[by.score$x < 4])))
  }
  wt.noise.library.genes<-as.vector(unique(temp.wt.emap$gene_name[abs(temp.wt.emap$score) > 2]))
  wt.noise.library.genes<-wt.noise.library.genes[!is.na(wt.noise.library.genes)]
  max.abs.score.temp.emap<-with(temp.emap, aggregate(score, by=list(gene_name=gene_name), abs.max))
  high.abs.score.subset<-as.vector(subset(max.abs.score.temp.emap, x>6)$gene_name)
  wt.noise.library.genes<-wt.noise.library.genes[! wt.noise.library.genes %in% high.abs.score.subset]
  low.score.library.genes<-unique(c(low.score.library.genes, wt.noise.library.genes))
  by.score.clean<-by.score[! by.score$Group.1 %in% low.score.library.genes,]
  by.score.clean.ordered<-by.score.clean[order(by.score.clean$x, decreasing=T),]
  by.score.ordered<-by.score[order(by.score$x, decreasing=T),]
  subunits<-as.vector(unique(by.score.clean.ordered$Group.1))
  final.emap.dim1<-data.frame()
  final.emap.dim2<-data.frame()
  for (j in 1:length(subunits)) {
    subunit = subunits[j]
    temp<-subset(temp.emap, temp.emap$gene_name==subunit)
    temp.dim1<-temp[match(ordered_genes.dim1, temp$Gene),]
    final.emap.dim1<-rbind(final.emap.dim1, temp.dim1)
    temp.dim2<-temp[match(ordered_genes.dim2, temp$Gene),]
    final.emap.dim2<-rbind(final.emap.dim2, temp.dim2)
  }
  for (k in 1:length(low.score.library.genes)) {
    low.subunit = low.score.library.genes[k]
    temp<-subset(temp.emap, temp.emap$gene_name==low.subunit)
    temp.dim1<-temp[match(ordered_genes.dim1, temp$Gene),]
    final.emap.dim1<-rbind(final.emap.dim1, temp.dim1)
    temp.dim2<-temp[match(ordered_genes.dim2, temp$Gene),]
    final.emap.dim2<-rbind(final.emap.dim2, temp.dim2)
  }
  final.emap.dim1<-cbind(final.emap.dim1, list("order"=seq(1, length(ordered_genes.dim1), 1)))
  final.emap.dim1<-cbind(final.emap.dim1, list("scaled.order"= multiscale.gene.dim1.order.table$dim1.order.scale))
  final.emap.dim2<-cbind(final.emap.dim2, list("order"=rev(seq(1,length(ordered_genes.dim2), 1))))
  final.emap.dim2<-cbind(final.emap.dim2, list("scaled.order"= multiscale.gene.dim2.order.table$dim2.order.scale))
  par(fig=c(0, 0.5, 0.35, 1)) #### MSD PLOT
  x<-fit$points[,1]
  y<-fit$points[,2]
  plot(x, y, type="n", xlab="", ylab="", main=compl, axes=F)
  box()
  text(x,y,labels=row.names(temp.emap.matrix), cex=0.6)
  par(fig=c(0,0.5,0,0.6), new=T)   #### TITRATION PLOT DIM 1
  plot(final.emap.dim1$scaled.order, final.emap.dim1$score, xaxt="n", type="n", xlab="", ylab="E-MAP score", ylim=c(emap.min,emap.max))
  if (length(na.omit(final.emap.dim1$score[final.emap.dim1$Gene == "GSP1-NAT" & final.emap.dim1$gene_name %in% subunits])) > 0) {
    abline(h=(min.na(final.emap.dim1$score[final.emap.dim1$Gene == "GSP1-NAT" & final.emap.dim1$gene_name %in% subunits])), lwd=0.5)
  } else {
    abline(h=0, lwd=0.5)
  }
  abline(v=final.emap.dim1$scaed.order[final.emap.dim1$Gene == "GSP1-NAT"], lwd=0.5)
  axis(1, at=multiscale.gene.dim1.order.table$dim1.order.scale, labels=ordered_genes.dim1, las=2, cex.axis=0.5)
  for (l in 1:length(low.score.library.genes)) {
    no.signal.subunit = low.score.library.genes[l]
    points(final.emap.dim1$scaled.order[final.emap.dim1$gene_name == no.signal.subunit], final.emap.dim1$score[final.emap.dim1$gene_name == no.signal.subunit], col = "grey", pch=1)
  }
  for (m in length(subunits):1) {
    subunit = subunits[m]
    points(final.emap.dim1$scaled.order[final.emap.dim1$gene_name == subunit], final.emap.dim1$score[final.emap.dim1$gene_name == subunit], pch=19, col=colors[m])
  }
  par(fig=c(0.4,0.9,0.35,1), new=T)   #### TITRATION PLOT DIM 2
  plot(final.emap.dim2$score, final.emap.dim2$scaled.order, yaxt="n", type="n", ylab="", xlab="E-MAP score", xlim=c(emap.min,emap.max))
  if (length(na.omit(final.emap.dim2$score[final.emap.dim2$Gene == "GSP1-NAT" & final.emap.dim2$gene_name %in% subunits])) > 0) {
    abline(v=(min.na(final.emap.dim2$score[final.emap.dim2$Gene == "GSP1-NAT" & final.emap.dim2$gene_name %in% subunits])), lwd=0.5)
  } else {
    abline(v=0, lwd=0.5)
  }
  abline(h=final.emap.dim2$scaled.order[final.emap.dim2$Gene == "GSP1-NAT"], lwd=0.5)
  axis(4, at=multiscale.gene.dim2.order.table$dim2.order.scale, labels=ordered_genes.dim2, las=2, cex.axis=0.5)
  for (l in 1:length(low.score.library.genes)) {
    no.signal.subunit = low.score.library.genes[l]
    points(final.emap.dim2$score[final.emap.dim2$gene_name == no.signal.subunit], final.emap.dim2$scaled.order[final.emap.dim2$gene_name == no.signal.subunit], col = "grey", pch=1)
  }
  for (m in length(subunits):1) {
    subunit = subunits[m]
    points(final.emap.dim2$score[final.emap.dim2$gene_name == subunit],final.emap.dim2$scaled.order[final.emap.dim2$gene_name == subunit], pch=19, col=colors[m])
  }
  par(fig=c(0.45,0.9,0,0.5),new=T)   #### LIBRAY GENE LEGEND
  plot(final.emap.dim1$scaled.order, final.emap.dim1$score,type = "n", axes = FALSE, ann = FALSE)
  if (length(subunits)>0) {
    if (length(subunits)<=15) {
      legend("topleft", pch=19, legend=subunits, col=colors[1:length(subunits)], cex=0.5)
    } else {
      legend("topleft", pch=19, legend=subunits[1:15], col=colors[1:length(subunits)][1:15], cex=0.5)
      legend("top", pch=19, legend=subunits[16:length(subunits)], col=colors[1:length(subunits)][16:length(subunits)], cex=0.5)
    }
  }
  library.genes.in.complexes<-append(library.genes.in.complexes, subunits)
  if (length(low.score.library.genes) > 0) {
    if (length(low.score.library.genes) <=15) {
      legend("topright", pch=1, legend=low.score.library.genes, col="grey", cex=0.5)
    } else {
      legend("bottom", pch=1, legend=low.score.library.genes[1:15], col="grey", cex=0.5)
      legend("topright", pch=1, legend=low.score.library.genes[16:length(low.score.library.genes)], col="grey", cex=0.5)
    }
  }
  library.genes.in.complexes<-append(library.genes.in.complexes, low.score.library.genes)
  ##### NOW REPEAT THE PLOTS WITHOUT SCALED ORDERING
  par(fig=c(0, 0.5, 0.35, 1)) #### MSD PLOT
  x<-fit$points[,1]
  y<-fit$points[,2]
  plot(x, y, type="n", xlab="", ylab="", main=compl, axes=F)
  box()
  text(x,y,labels=row.names(temp.emap.matrix), cex=0.6)
  par(fig=c(0,0.5,0,0.6), new=T)   #### TITRATION PLOT DIM 1
  plot(final.emap.dim1$order, final.emap.dim1$score, xaxt="n", type="n", xlab="", ylab="E-MAP score", ylim=c(emap.min,emap.max))
  if (length(na.omit(final.emap.dim1$score[final.emap.dim1$Gene == "GSP1-NAT" & final.emap.dim1$gene_name %in% subunits])) > 0) {
    abline(h=(min.na(final.emap.dim1$score[final.emap.dim1$Gene == "GSP1-NAT" & final.emap.dim1$gene_name %in% subunits])), lwd=0.5)
  } else {
    abline(h=0, lwd=0.5)
  }
  abline(v=final.emap.dim1$x[final.emap.dim1$Gene == "GSP1-NAT"], lwd=0.5)
  axis(1, at=1:length(ordered_genes.dim1), labels=ordered_genes.dim1, las=2, cex.axis=0.5)
  for (l in 1:length(low.score.library.genes)) {
    no.signal.subunit = low.score.library.genes[l]
    points(final.emap.dim1$order[final.emap.dim1$gene_name == no.signal.subunit], final.emap.dim1$score[final.emap.dim1$gene_name == no.signal.subunit], col = "grey", pch=1)
  }
  for (m in length(subunits):1) {
    subunit = subunits[m]
    points(final.emap.dim1$order[final.emap.dim1$gene_name == subunit], final.emap.dim1$score[final.emap.dim1$gene_name == subunit], pch=19, col=colors[m])
  }
  par(fig=c(0.4,0.9,0.35,1), new=T)   #### TITRATION PLOT DIM 2
  plot(final.emap.dim2$score, final.emap.dim2$order, yaxt="n", type="n", ylab="", xlab="E-MAP score", xlim=c(emap.min,emap.max))
  if (length(na.omit(final.emap.dim2$score[final.emap.dim2$Gene == "GSP1-NAT" & final.emap.dim2$gene_name %in% subunits])) > 0) {
    abline(v=(min.na(final.emap.dim2$score[final.emap.dim2$Gene == "GSP1-NAT" & final.emap.dim2$gene_name %in% subunits])), lwd=0.5)
  } else {
    abline(v=0, lwd=0.5)
  }
  abline(h=final.emap.dim2$order[final.emap.dim2$Gene == "GSP1-NAT"], lwd=0.5)
  axis(4, at=length(ordered_genes.dim2):1, labels=ordered_genes.dim2, las=2, cex.axis=0.5)
  for (l in 1:length(low.score.library.genes)) {
    no.signal.subunit = low.score.library.genes[l]
    points(final.emap.dim2$score[final.emap.dim2$gene_name == no.signal.subunit], final.emap.dim2$order[final.emap.dim2$gene_name == no.signal.subunit], col = "grey", pch=1)
  }
  for (m in length(subunits):1) {
    subunit = subunits[m]
    points(final.emap.dim2$score[final.emap.dim2$gene_name == subunit],final.emap.dim2$order[final.emap.dim2$gene_name == subunit], pch=19, col=colors[m])
  }
  par(fig=c(0.45,0.9,0,0.5),new=T)   #### LIBRAY GENE LEGEND
  plot(final.emap.dim1$scaled.order, final.emap.dim1$score,type = "n", axes = FALSE, ann = FALSE)
  if (length(subunits)>0) {
    if (length(subunits)<=15) {
      legend("topleft", pch=19, legend=subunits, col=colors[1:length(subunits)], cex=0.5)
    } else {
      legend("topleft", pch=19, legend=subunits[1:15], col=colors[1:length(subunits)][1:15], cex=0.5)
      legend("top", pch=19, legend=subunits[16:length(subunits)], col=colors[1:length(subunits)][16:length(subunits)], cex=0.5)
    }
  }
  if (length(low.score.library.genes) > 0) {
    if (length(low.score.library.genes) <=15) {
      legend("topright", pch=1, legend=low.score.library.genes, col="grey", cex=0.5)
    } else {
      legend("bottom", pch=1, legend=low.score.library.genes[1:15], col="grey", cex=0.5)
      legend("topright", pch=1, legend=low.score.library.genes[16:length(low.score.library.genes)], col="grey", cex=0.5)
    }
  }
}
dev.off()

