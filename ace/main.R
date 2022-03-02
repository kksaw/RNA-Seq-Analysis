# only the main commands - for detailed description, check individual algorithm files

setwd("C:/Users/khaik/Documents/camb/sawarkar lab/ace")

# ----- SPECIAL FUNCTIONS -----
which(sapply(gene, function(x) any(grepl("ENSMUSG00000056458", x))))  # return counts in each sample for a specific gene)

assayNames(gse)  # view name of assays (ie. counts, abundance, length)
head(assay(gse),3)  # assay(gse) return counts
colSums(assay(gse))
rowRanges(gse)  # shows ranges for first 5 and last 5 genes
seqinfo(rowRanges(gse))  # metadata abt sequences (ie. chromosones, seqlength, isCircular, genome)
colData(gse) # info from sample_table.csv

# results filter
# many genes w DE at false discovery rate (FDR) of 10%, so
# - lower FDR TH (padj)
# - increase L2FC TH (lfcThreshold) from 0
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold = 1) # needs to be more than 2x increase/decrease
table(resLFC1$padj < 0.1)

#! compare between two celltypes
# results(dds, contrast = c("cell", "Wt", "H14"))

# multiple testing - BH adjusted p value
#! if we don't use adjusted p values, there are 9536/33330 genes with p value < 0.05 which is a lot higher than expected 5% = 1666
#! ie. 1666/9536 = 17.4% false positives
sum(res$pvalue < 0.05, na.rm=TRUE)  #! 9536
sum(!is.na(res$pvalue)) #! 33330

# if we say 10% false positives is acceptable
sum(res$padj < 0.1, na.rm=TRUE)  # 8607

# subset this then sort by l2fc to get significant genes with strongest down regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ]) # strongest downregulation
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])  # strongest upregulation

# interactive MDS plot
labels <- paste(coldata$names, type, status)
group <- paste(type, status, sep=".")
group <- factor(group)
glMDSPlot(dgeObj, labels=labels, groups=group, folder="mds")



# START HERE
# ----- LOAD DATA -----
library("tidyverse")
library("AnnotationDbi")
library("org.Hs.eg.db")
csvfile <- "sample details/sample_info time series.csv"
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = TRUE)
coldata$files <- str_c("quants_time series/", rownames(coldata), "/quant.sf")
file.exists(coldata$files)

sample_info <- coldata[is.element(coldata$Sample, c(438,439,440)), ]
sample_info <- sample_info[!is.element(sample_info$Day, 14), ]

genefile <- "proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv"
genedata <- read.csv(genefile, row.names = 1, stringsAsFactors = TRUE)
genedata <- genedata[, c(1,3,5,15)]
levels(genedata$Pathway) <- cbind(levels(genedata$Pathway), 'NA')
mypalette = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
              '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', 
              '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffffff')
pcol = cbind(levels(genedata$Pathway), mypalette)
rownames(pcol) <- pcol[,1]
pcol <- pcol[,2]

library("tximeta")
library("MatrixGenerics")
library("SummarizedExperiment")
se <- tximeta(sample_info)
gse <- summarizeToGene(se)
levels(gse$Sample <- factor(gse$Sample))
levels(gse$Day <- factor(gse$Day))

# ----- FILTERING -----
# check if duplicates and remove (none)
o <- order(rowSums(assay(gse)), decreasing=TRUE)
gse <- gse[o,]
d <- duplicated(rownames(gse))
gse <- gse[!d]
nrow(gse)
rm(o); rm(d)

# remove genes with low count (choose one)
# 1. remove genes with 0/1 count across all samples
keep <-  rowSums(assay(gse)) > 1
# 2. remove genes if not at least 3 samples with count of 10 or higher
# keep <- rowSums(counts(dds) >= 10) >=3
# 3. by cpm
# myCPM <- cpm(assay(gse))
# plot(myCPM[,1], assay(gse)[,1], ylim=c(0,50), xlim=c(0,3)) # choose a cpm that corresponds to 10-15 counts
# abline(v=0.25) # threshold = 0.25
# thresh <- myCPM > 0.25
# table(rowSums(thresh)) # summary of how many trues in each row (here 12380 genes are < 0.5 in all 12 samples)
# keep <- rowSums(thresh) >= 3 # 3 replicates
# # summary(keep)

gse <- gse[keep,]
rownames(gse) <- sub("\\..*","",rownames(gse))
ens.str <- rownames(gse)
nrow(gse)

library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)

# ----- DESEQ2 -----
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ Sample + Day)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)  # compares last variable in design formula
res$symbol <- mapIds(org.Hs.eg.db, keys=ens.str,
                     column="SYMBOL", keytype="ENSEMBL", multiVals="first")
restable <- as.data.frame(res[order(res$padj),])

res4 <- results(dds, name="Day_4_vs_1")
res8 <- results(dds, name="Day_8_vs_1")
#res14 <- results(dds, name="Day_14_vs_1")
res28 <- results(dds, name="Day_28_vs_1")
res63 <- results(dds, name="Day_63_vs_1")
res84 <- results(dds, name="Day_84_vs_1")

a <-  gse[is.element(rownames(gse), genedata$Human_gene_ID),]

res4$symbol <- mapIds(org.Hs.eg.db, keys=ens.str,
                      column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res4table <- as.data.frame(res4[order(res4$padj),])

# ----- MA PLOT -----
make_MAplot <- function(day, res, yl=5) {
  resLFC <- lfcShrink(dds, coef=paste("Day_", day, "_vs_1", sep=""), type="apeglm")
  plotMA(resLFC, ylim = c(-yl,yl))
}


# ---- COUNT PLOT -----
dds <- estimateSizeFactors(dds)
normalizationFactors(dds)
#a <- counts(dds, normalized = TRUE)

day = "4"
res = res4

rownames(dds) <- mapIds(org.Hs.eg.db, keys=rownames(dds),
                        column="SYMBOL", keytype="ENSEMBL", multiVals="first")
dds$Day <- as.numeric(as.character(dds$Day))

make_countplot <- function(ddsc, thisgene) {
  ggplot(ddsc, aes(x = Day, y = count, color = Sample, group = Sample)) +
    geom_point() + stat_summary(fun = mean, geom = "line") +
    scale_y_log10() + labs(title = thisgene, x="", y="") + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
}

library("cowplot")
make_countplot_combined <- function(day, res, ngenes=12) {
  par(mfrow=c(3,4))
  
  # all genes
  topGenes <- head(order(res$padj), ngenes)
  plot_list <- list()
  
  for (i in topGenes) {
    thisgene = rownames(dds)[i]
    ddsc <- plotCounts(dds, i, intgroup = c("Day", "Sample"), returnData=TRUE)
    plot_list <- append(plot_list, list(make_countplot(ddsc, thisgene)))
  }
  
  plot_grid(plotlist = plot_list)
  ggsave(paste("fig/Day",day,"_top",toString(ngenes),"_all.png", sep=""), width=10, height=5)
  
  # only PTSS genes
  res_sel <- res[is.element(rownames(dds), genedata$Human_Symbol), ]
  topGenes_pos <- head(order(res_sel$padj), ngenes)
  topGenes <- mapIds(org.Hs.eg.db, keys=rownames(res_sel),
                     column="SYMBOL", keytype="ENSEMBL", multiVals="first")[topGenes_pos]
  plot_list <- list()
  
  for (i in topGenes) {
    ddsc <- plotCounts(dds, i, intgroup = c("Day", "Sample"), returnData=TRUE)
    plot_list <- append(plot_list, list(make_countplot(ddsc, i)))
  }
  
  plot_grid(plotlist = plot_list)
  ggsave(paste("fig/Day",day,"_top",toString(ngenes),"_ptss.png", sep=""), width=10, height=5)
}

make_countplot_combined("8", res8)
make_countplot_combined("28", res28)
make_countplot_combined("63", res63)
make_countplot_combined("84", res84)

# ----- HEATMAP ------
library("pheatmap")
library("grid")
betas <- coef(dds)
colnames(betas)

#gpathway <- data.frame(pathway = genedata[head(order(match(genedata$Human_Symbol, rownames(mat))),50), ]$Pathway)

day = "84"
res = res84

make_gpathway <- function(glist) {
  gpathway <- plyr::join(data.frame(Human_Symbol = glist), genedata[,c(1,2)], type="left")
  gpathway[is.na(gpathway)] <- "NA"
  gpathway <- data.frame(pathway = gpathway[,2])
  rownames(gpathway) <- glist
  return(gpathway)
}

plot_heatmap <- function(topGenes, thr, day, method) {
  mat <- betas[topGenes, -c(1,2,3)]
  rownames(mat) <- mapIds(org.Hs.eg.db, keys=rownames(mat),
                          column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  
  mat[mat < -thr] <- -thr
  mat[mat > thr] <- thr
  
  gpathway <- make_gpathway(rownames(mat))
  my_color = list(pathway = pcol[unique(gpathway$pathway)])
  
  pheatmap(mat, cluster_col = FALSE, cluster_rows = FALSE, 
           annotation_row = gpathway, annotation_colors = my_color,
           cellwidth= 12, cellheight = 8, fontsize = 8, filename=paste("fig/Day",day,"vs1_byPadj_",method,".png", sep=""))
  
  hmap <- pheatmap(mat, cluster_col=FALSE, cluster_rows = TRUE)
  gpathway <- make_gpathway(hmap$gtable$grobs[[4]]$label)
  my_color = list(pathway = pcol[unique(gpathway$pathway)])
  
  pheatmap(mat, cluster_col = FALSE, cluster_rows = TRUE, 
           annotation_row = gpathway, annotation_colors = my_color,
           cellwidth= 12, cellheight = 8, fontsize = 8, filename=paste("fig/Day",day,"vs1_",method,".png", sep="")) 
}

heatmap_fn <- function(day, res, thr1=7, thr2=5) {
  #including non PTSS genes
  topGenes <- head(order(res$padj), 50)
  plot_heatmap(topGenes, thr1, day, method="all")
  
  #only PTSS genes
  gnames <- mapIds(org.Hs.eg.db, keys=rownames(res),
                   column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  res_sel <- res[is.element(gnames, genedata$Human_Symbol), ]
  topGenes_pos <- head(order(res_sel$padj), 50)
  topGenes <- rownames(res_sel)[topGenes_pos]
  plot_heatmap(topGenes, thr2, day, method="ptss")
}

# ----- GRIDPLOT (not done) -----
library(gplots)
plot_grid <- function() {
  gnames <- mapIds(org.Hs.eg.db, keys=rownames(res[order(res$padj),]),
              column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  res_sel <- res[is.element(gnames, genedata$Human_Symbol), ]
  is_ptss <- plyr::join(data.frame(Human_Symbol = gnames), genedata[,c(1,2)], type="left")
  is_ptss[is.na(is_ptss)] <- "NA"
  
  pno <- cbind(levels(genedata$Pathway), 1:18)
  rownames(pno) <- pno[,1]
  pno <- pno[,2]
  a <- plyr::revalue(is_ptss$Pathway, pno)
  a1 <- as.data.frame(apply(array(head(a,1000), dim=c(20,50)), 2, as.numeric))
  a1 <- t(a1)
  
  p <- pheatmap(a1, cluster_col = FALSE, cluster_row = FALSE, color = mypalette, filename="test.png")
  N_col_areas <- length(p$gtable$grobs[[3]]$children[[1]]$y)
  N_labels <- length(p$gtable$grobs[[3]]$children[[2]]$label)
  yseq1 <- seq(0,1,length.out=N_col_areas)
  center <- yseq1[2]/2
  yseq2 <- seq(center, 1-center, length.out=N_labels)
  
  for (k in 1:(N_labels)) {
    p$gtable$grobs[[3]]$children[[2]]$y[[k]] <- 
      sum(yseq2[k]*min(unit(1, "npc"), unit(150, "bigpts")),
          unit(1, "npc"), -1*min(unit(1, "npc"), unit(150, "bigpts")))
  }  
  plot.new()
  png("test.png")
  p
  dev.off()
  
  gpathway <- make_gpathway(hmap$gtable$grobs[[4]]$label)
  my_color = list(pathway = pcol[unique(gpathway$pathway)])
  
  pheatmap(a1, cluster_col = FALSE, cluster_rows = FALSE, 
           annotation_colors = my_color,
           cellwidth= 12, cellheight = 8, fontsize = 8, filename=paste("fig/Day",day,"vs1_",method,".png", sep="")) 
  library(gplots)
  heatmap.2(a1, breaks=1:18, col=mypalette)
  
  

}

ggplot(x, 
       aes(x = Name, y = variable, fill = factor(value))) + 
  geom_tile() + 
  scale_fill_manual(values=colors)+
  scale_x_discrete(name="Time Period", limits= rownames(df))

# ----- SAVECSV -----
save_csv <- function(day, res) {
  gnames <- mapIds(org.Hs.eg.db, keys=rownames(res),
                   column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  res_sel <- res[is.element(gnames, genedata$Human_Symbol), ]
  rownames(res_sel) <- mapIds(org.Hs.eg.db, keys=rownames(res_sel),
                       column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  res_sel_table <- as.data.frame(res_sel[order(res_sel$padj), c(1,2,6)])
  write.csv(res_sel_table, file = paste("results/Day",day,"vs1.csv", sep=""))
}


# res <- results(dds, contrast=c("lps","trt","untrt"))  # to make sure we are comparing different treatments

summary(res)
res$symbol <- mapIds(org.Hs.eg.db, keys=ens.str,
                     column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, keys=ens.str,
                     column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res$genename <- mapIds(org.Hs.eg.db, keys=ens.str,
                       column="GENENAME", keytype="ENSEMBL", multiVals="first")

restable <- as.data.frame(res[order(res$padj),])
write.csv(results_deseq2, file = "results_deseq2.csv")


vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
dds <- estimateSizeFactors(dds)

#check normalisation
normalizationFactors(dds)

library("ggplot2")
# VISUALISATION
# Euclidean and Poisson distance
sampleDists <- dist(t(assay(vsd))) # Euclidean distance

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$lps, vsd$expt, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$lps, dds$expt, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA
plotPCA(vsd, intgroup = c("lps", "expt"))
pcaData <- plotPCA(vsd, intgroup = c( "cell", "expt"), returnData = TRUE)
pcaData

# nicer plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = cell, shape = expt)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


# Generalised PCA
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$cell <- dds$cell
gpca.dat$lps <- dds$lps
gpca.dat$expt <- dds$expt

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = cell, shape = expt)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# MDS
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = lps, shape = expt)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = lps, shape = expt)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")


# ----- edgeR -----
library(edgeR)

dgeObj <- DGEList(assay(gse))
barplot(dgeObj$samples$lib.size, names=paste(sample_info$Sample, "_Day", sample_info$Day, sep=""), las=2) #check if any discrepancies btw samples
title("Barplot of library sizes")

dgeObj$genes$symbol <- mapIds(org.Mm.eg.db, keys=ens.str,
                              column="SYMBOL", keytype="ENSEMBL", multiVals="first")
dgeObj$genes$entrez <- mapIds(org.Mm.eg.db, keys=ens.str,
                              column="ENTREZID", keytype="ENSEMBL", multiVals="first")
dgeObj$genes$genename <- mapIds(org.Mm.eg.db, keys=ens.str,
                              column="GENENAME", keytype="ENSEMBL", multiVals="first")
dgeObj$genes <- dgeObj$genes[, -c(2:7)] # delete col

dgeObj <- calcNormFactors(dgeObj)
# design <- model.matrix(~ pairing + condition) # for coldata_Wt/H14
design <- model.matrix(~ pairing + celltype) # for coldata_untrt/trt

dgeObj <- estimateDisp(dgeObj)
dgeObj$common.dispersion
plotBCV(dgeObj)

fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit)
topTags(lrt)

summary(de <- decideTests(lrt))
plotMD(lrt)
results_edgeR <- as.data.frame(topTags(lrt, n=Inf))
write.csv(results_edgeR, file="results_edgeR.csv", row.names=TRUE)

# ----- limma-voom -----
library(limma)
library(Glimma)
v <- voom(dgeObj, design, plot = TRUE)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")

# limma - fit linear model
fit_limma <- lmFit(v, design)
fit.cont <- contrasts.fit(fit_limma, coef=4)
fit.cont <- eBayes(fit.cont)
plotSA(fit.cont, main="Final model: Mean-variance trend")
summary(summa.fit <- decideTests(fit.cont))
results_limma <- as.data.frame(topTreat(fit.cont, n=Inf))
write.csv(results_limma, file="results_limma.csv", row.names=TRUE)


# check similarity 
setwd("C:/Users/khaik/Documents/camb/sawarkar lab/entrega/results/")
results_deseq2 <- read.csv("results_deseq2.csv", row.names = 1, stringsAsFactors = FALSE)
results_edgeR <- read.csv("results_edgeR.csv", row.names = 1, stringsAsFactors = FALSE)
results_limma <- read.csv("results_limma.csv", row.names = 1, stringsAsFactors = FALSE)

gdeseq2 <- rownames(results_deseq2)
gedgeR <- rownames(results_edgeR)
glimma <- rownames(results_limma)

n <- 1000
int <- Reduce(intersect, list(gdeseq2[1:n],gedgeR[1:n],glimma[1:n]))

gdetails <- results_edgeR[is.element(gedgeR, int), c(1,3,2,4,6)]
sdeseq2 <- results_deseq2[is.element(gdeseq2, int), c(2,5,6)]
sedgeR <- results_edgeR[is.element(gedgeR, int), c(5,7,8,9)]
slimma <- results_limma[is.element(glimma, int), c(5,8,9,10)]
results_combined_Wt_H14 <- cbind(gdetails, sdeseq2, sedgeR, slimma)
write.csv(results_combined,file="results_H14_only.csv",row.names=TRUE)

n <- 100
dboth <- intersect(rownames(results_combined_Wt)[1:n], rownames(results_combined_H14)[1:n])
dWt <- setdiff(rownames(results_combined_Wt)[1:n], rownames(results_combined_H14)[1:n])
dH14 <- setdiff(rownames(results_combined_H14)[1:n], rownames(results_combined_Wt)[1:n])


# volcano plot
library("RColorBrewer")
library("ggrepel")
make_volcanoplot <- function(day, res, ylim10=5, ylim1=80, ylim2=60, xlim1=10, xlim2=5, nlabel=12, n=1000) {
  restable <- as.data.frame(res[order(res$padj),c(1,2,6)])
  gnames <- mapIds(org.Hs.eg.db, keys=rownames(restable), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  mat <- restable[1:n,]
  
  mat$gene_family[is.element(gnames[1:n], genedata$Human_Symbol)] <- "proteostasis"
  mat$gene_family[!is.element(gnames[1:n], genedata$Human_Symbol)] <- "others"
  colnames(mat)[2] <- "L2FC"
  
  mat$topGenes <- ""
  mat$topGenes[1:nlabel] <- gnames[1:nlabel]
  ggplot(mat, aes(x = L2FC, y=-log10(padj), label=topGenes)) +
    geom_point(aes(colour=gene_family, size=log2(baseMean))) +
    ggtitle(paste("Day", day, "vs1", sep="")) + xlab("L2FC") + ylab("-log10 adjusted p-value") + 
    ylim(ylim10,ylim1) + xlim(-xlim1,xlim1) +
    scale_color_brewer(palette="Pastel1") + geom_text_repel()  
  
  ggsave(paste("fig/Day",day,"_volcano_all.png", sep=""), width=10, height=5)
  
  mat <- restable[is.element(gnames, genedata$Human_Symbol),][1:n,]
  colnames(mat)[2] <- "L2FC"
  mat$topGenes <- ""
  mat$topGenes[1:nlabel] <- gnames[is.element(gnames, genedata$Human_Symbol)][1:nlabel]
  ggplot(mat, aes(x = L2FC, y=-log10(padj), label=topGenes)) +
    geom_point(colour="#B3CDE3", aes(size=log2(baseMean))) +
    ggtitle(paste("Day", day, "vs1", sep="")) + xlab("L2FC") + ylab("-log10 adjusted p-value") + 
    ylim(0,ylim2) + xlim(-xlim2,xlim2) + geom_text_repel()  
  
  ggsave(paste("fig/Day",day,"_volcano_ptss.png", sep=""), width=10, height=5)
  
}



# ----- GENE CLUSTERING -----
library("genefilter")
mat <- assay(vsd[is.element(rownames(vsd), int)])
mat <- mat - rowMeans(mat) # see how much a gene deviates from gene's average across samples rather than absolute expression strength
anno <- as.data.frame(sample_info[, c("lps", "expt")])
mat.str <- sub("\\..*","",rownames(mat))
rownames(mat) <- mapIds(org.Mm.eg.db, keys=mat.str, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
pheatmap(mat, annotation_col = anno, cluster_rows = F)



# ----- PLOTTING -----
# 1. counts plot for individual gene
topGene <- rownames(res)[which.min(res$padj)] # lowest FDR TH
plotCounts(dds, gene = topGene, intgroup=c("lps"))

# nicer plot - normalised counts for a single gene over treatment group
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("lps", "expt"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = lps, y = count, color = expt)) + 
  scale_y_log10() + geom_beeswarm(cex = 3)

# with a line joining same pairing - something wrong here
ggplot(geneCounts, aes(x = lps, y = count, color = expt, group = expt)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# 2. MA-plot aka mean-difference plot/Bland-Altman plot
library("apeglm")
resultsNames(dds)
# shrink l2fc using apeglm (good for noisy LFC estimates while giving low bias LFC estimates for true large differences)
res <- lfcShrink(dds, coef="cell_H14_vs_Wt", type="apeglm") 

# MA plot: genes with p value < 0.1 (default) will be shown in red
DESeq2::plotMA(res, ylim = c(-10,10))

# label individual points - circle and label interested gene
plotMA(res, ylim = c(-10,10))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# histogram of p values - exclude genes with very small counts which would otherwise generate spikes in histogram
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# edgeR
detags <- rownames(dgeObj)[as.logical(de)]
par(mfrow=c(1,2))
plotSmear(lrt.TvsU, de.tags=detags)
signif <- -log10(results_edgeR$FDR)
plot(results_edgeR$logFC, signif, pch=16)
points(results_edgeR[detags,"logFC"],-log10(results_edgeR[detags,"FDR"]),pch=16,col="red")
