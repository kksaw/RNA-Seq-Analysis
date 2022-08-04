# only the main commands - for detailed description, check individual algorithm files

setwd("C:/Users/khaik/Documents/camb/sawarkar lab/entrega")

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
library("AnnotationDbi")
library("org.Mm.eg.db")
csvfile <- "sample_info.csv"
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
file.exists(coldata$files)  # check if quant files exist

coldata_Wt <- coldata[coldata$cell == "Wt", ]
coldata_H14 <- coldata[coldata$cell == "H14", ]
coldata_untrt <- coldata[coldata$lps == "untrt",]
coldata_trt <- coldata[coldata$lps == "trt", ]

sample_info <- coldata_untrt

library("tximeta")
library("MatrixGenerics")
library("SummarizedExperiment")
se <- tximeta(sample_info)  # tximeta expects at least 2 columns in coldata - files (pointer to quant.sf files), names (to identify sample); download annotation data at transcript level
gse <- summarizeToGene(se) # summarise transcript-level quantifications to gene-level

levels(gse$cell <- factor(gse$cell))
levels(gse$lps <- factor(gse$lps))
levels(gse$expt <- factor(gse$expt))

# R prefers first level of factor to be ref (ie. control, wild type, untreated) sample. If colData table was assembled with treated samples set as reference, we need to use relevel like below. If not, skip this step
library("magrittr")
celltype <- gse$cell %<>% relevel("Wt")  # dont run if coldata_Wt/H14
pairing <- gse$expt
condition <- gse$lps %<>% relevel("untrt") # dont run if coldata_trt/untrt

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
nrow(gse)

ens.str <- sub("\\..*","",rownames(gse))

# ----- DESEQ2 -----
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ expt + cell)  #! 54309 genes
vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)

library("dplyr")
library("ggplot2")
dds <- estimateSizeFactors(dds)

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


# DIFFERENTIAL EXPRESSION
dds <- DESeq(dds)
results_deseq2 <- results(dds)  # compares last variable in design formula
# res <- results(dds, contrast=c("lps","trt","untrt"))  # to make sure we are comparing different treatments

summary(results_deseq2)
results_deseq2$symbol <- mapIds(org.Mm.eg.db, keys=ens.str,
                     column="SYMBOL", keytype="ENSEMBL", multiVals="first")
results_deseq2$entrez <- mapIds(org.Mm.eg.db, keys=ens.str,
                     column="ENTREZID", keytype="ENSEMBL", multiVals="first")
results_deseq2$genename <- mapIds(org.Mm.eg.db, keys=ens.str,
                       column="GENENAME", keytype="ENSEMBL", multiVals="first")

results_deseq2 <- as.data.frame(results_deseq2[order(results_deseq2$padj),])
write.csv(results_deseq2, file = "results_deseq2.csv")

# ----- edgeR -----
library(edgeR)

dgeObj <- DGEList(gse)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2) #check if any discrepancies btw samples
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
mat_wt <- results_combined_Wt[1:n, c(3,5,6,8)]
mat_wt$cell[is.element(rownames(mat_wt), dWt)] <- "wt"
mat_wt$cell[is.element(rownames(mat_wt), dboth)] <- "both"
mat_h14 <- results_combined_H14[1:n, c(3,5,6,8)]
mat_h14$cell[is.element(rownames(mat_h14), dH14)] <- "h14"
mat_h14$cell[is.element(rownames(mat_h14), dboth)] <- "both"
mat <- rbind(mat_wt, mat_h14)
mat <- mat[order(mat$padj),]
colnames(mat)[3] <- "logFC"
mat$padj[1:3] <- 1e-320

library("ggrepel")
mat$top20 <- ""
mat$top20[1:15] <- mat$symbol[1:15]
ggplot(mat, aes(x = logFC, y=-log10(padj), label=top20)) +
  geom_point(aes(colour=cell, size=logCPM)) +
  ggtitle("DEGs") + xlab("logFC") + ylab("-log10 adjusted p-value") + ylim(10,330) +
  scale_color_brewer(palette="Pastel1") + geom_text_repel()

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
