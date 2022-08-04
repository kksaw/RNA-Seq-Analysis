setwd("C:/Users/khaik/Documents/camb/sawarkar lab/entrega")

library(ggplot2)

# ----- functions -----
load_count_data <- function(sample_info, filename) {
  Rex <- read.delim("quants_eisa/ExonicCounts.txt")
  Rin <- read.delim("quants_eisa/IntronicCounts.txt")
  
  sample_name <- gsub(" ", "", sample_info$sample, fixed = TRUE)
  
  if (filename == "2 H14_trt vs untrt") {
    ref <- sample_name[4:6]
    comp <- sample_name[1:3]
  } else {
    ref <- sample_name[1:3]
    comp <- sample_name[4:6]
  }
  
  if (is.element(filename, c("1 Wt_trt vs untrt", "2 H14_trt vs untrt"))){
    cond <- sample_info$lps
  } else {
    cond <- sample_info$cell    
  }
  
  Rex <- Rex[, sample_name]
  Rin <- Rin[, sample_name]
  
  return(list("Rex" = Rex, "Rin" = Rin, "sample_name" = sample_name, 
              "cond" = cond, "ref" = ref, "comp" = comp))
}

prepare_y <- function (Rex, Rin, comp, ref, th=3.0) {
  # 1. Normalisation
  Rall <- Rex + Rin
  fracIn <- colSums(Rin)/colSums(Rall)
  print(summary(fracIn))
  
  # scale counts to the mean library size, separately for exons and introns
  Nex <- t(t(Rex) / colSums(Rex) * mean(colSums(Rex)))
  Nin <- t(t(Rin) / colSums(Rin) * mean(colSums(Rin)))
  
  # log transform (add a pseudocount of 8)
  NLex <- log2(Nex + 8)
  NLin <- log2(Nin + 8)
  
  # 2. Identify quantifiable genes
  quantGenes <- rownames(Rex)[ rowMeans(NLex) > th & rowMeans(NLin) > th ]
  print(paste(length(quantGenes), " quantifiable genes", sep=""))
  
  # 3. Calculate changes in introns ((Din) and exons (Dex) and differece (Dex - Din)
  Dex <- NLex[, comp] - NLex[, ref]
  Din <- NLin[, comp] - NLin[, ref]
  Dex.Din <- Dex - Din 
  
  return(list("fracIn" = fracIn, "quantGenes" = quantGenes, 
              "Dex" = Dex, "Din" = Din, "Dex.Din" = Dex.Din))
}

DE_analysis <- function(cnt, quantGenes, filename, 
                        cnt_type = "EI", save_flag = FALSE) {
  
  y <- DGEList(counts = cnt, genes = data.frame(ENTREZID = rownames(cnt)))
  
  # select quantifiable genes and normalize
  y <- y[quantGenes, ]
  y <- calcNormFactors(y)
  
  y$genes$symbol <- mapIds(org.Mm.eg.db, keys=rownames(y),
                           column="SYMBOL", keytype="ENTREZID", multiVals="first")
  y$genes$genename <- mapIds(org.Mm.eg.db, keys=rownames(y),
                             column="GENENAME", keytype="ENTREZID", multiVals="first")
  
  # design matrix with interaction term
  if (substr(cnt_type,1,2) == "EI") {
    region <- factor(c("ex","ex","ex","ex","ex","ex","in","in","in","in","in","in"), levels = c("in", "ex"))
    cond <- rep(cond, 2)
    design <- model.matrix(~ region * cond)
    rownames(design) <- colnames(cnt)
    savefile = paste("/results_edgeR_", cnt_type, ".csv", sep = "")
  } else {
    design <- model.matrix(~ cond)
    rownames(design) <- colnames(cnt)
    if (cnt_type == "exon") {
      savefile = "/results_edgeR_exon.csv"
    } else if (cnt_type == "intron") {
      savefile = "/results_edgeR_intron.csv"
    }
  }
  
  # estimate model parameters and calculate LLR btw full and reduced models
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  
  # results table
  res_edgeR <- as.data.frame(topTags(lrt, n=Inf))
  print(y$samples$lib.size)
  
  if (save_flag) {
    write.csv(res_edgeR, file=paste("results_eisa/", filename, savefile, sep=""), row.names=FALSE)
  }
  
  return(list("y" = y, "res_edgeR" = res_edgeR))  
}

plot_libsize <- function(y, cond, ref, comp, type="exon"){
  libsize <-data.frame(Experiment = c("A", "B", "C"), 
                       Size = y$samples$lib.size,
                       Condition = cond)
  
  ggplot(data = libsize, aes(x = Experiment, y = Size, fill = Condition)) +
    geom_bar(stat="identity", position = position_dodge()) +
    labs(title = paste("Barplot of library sizes - ", type, sep="")) +
    scale_fill_brewer(palette = "Blues") +
    theme_classic()
  
  ggsave(paste("results_eisa/", filename, "/libsize_", type, ".png", sep=""), width=5, height=3)
}

print_summary <- function(res_df, ori_df, quantGenes, type="edgeR_exon"){
  if (substr(type,1,5) == "edgeR"){
    print(paste("Total genes included: ", nrow(res_df), "   ",
                "Up: ", sum(res_df$logFC > 0), "   ",
                "Down: ", sum(res_df$logFC < 0), "   ",
                "No change: ", sum(res_df$logFC == 0), 
                sep = ""))
    print(paste("Mean change: ", mean(ori_df[quantGenes,])))
    print(paste("Correlation: ", cor(ori_df[quantGenes,1], ori_df[quantGenes,2]), cor(ori_df[quantGenes,1], ori_df[quantGenes,3]), cor(ori_df[quantGenes,2], ori_df[quantGenes,3])))
    
    print("FDR < 0.05")
    a <- res_df[res_df$FDR < 0.05, ]
    print(paste("Total genes included: ", nrow(a), "   ",
                "Up: ", sum(a$logFC > 0), "   ",
                "Down: ", sum(a$logFC < 0), "   ",
                "No change: ", sum(a$logFC == 0), 
                sep = ""))
    a1 <- rownames(a)
    a2 <- ori_df[is.element(rownames(ori_df), a1),]
    print(paste("Mean change: ", mean(a2)))
    print(paste("Correlation: ", cor(a2[,1], a2[,2]), cor(a2[,1], a2[,3]), cor(a2[,2], a2[,3])))
    
    
    if (substr(type, 7, 8) == "EI"){
      g <- c("<0", ">0")
      c <- c("#BCBDDC", "#FDAE6B", "#FFED6F")
      
    } else {
      g <- c("<0", "0", ">0")
      f <- c("All", "All", "All", "< 0.05", "< 0.05", "< 0.05")
      df <-data.frame(logFC = factor(g, levels = unique(g)),
                      ngenes = c(sum(res_df$logFC < 0) - sum(a$logFC < 0), sum(res_df$logFC == 0) - sum(a$logFC == 0), sum(res_df$logFC > 0) - sum(a$logFC > 0),
                                 sum(a$logFC < 0), sum(a$logFC == 0), sum(a$logFC > 0)),
                      FDR = factor(f, levels = unique(f)))
      
      ggplot(data = df, aes(x = logFC, y = ngenes, fill = logFC)) +
        geom_bar(stat="identity", position = "stack", aes(alpha = FDR)) +
        labs(title = paste("Delta", substr(type, 6, nchar(type)), sep="")) +
        ylim(0,9000) +
        scale_fill_manual(values = c("#9ECAE1", "#A1D99B", "#FC9272")) +
        scale_alpha_manual(values = c(0.95, 0.5)) +
        theme_classic()
      
      ggsave(paste("results_eisa/", filename, "/libsize_", type, ".png", sep=""), width=3, height=3)
    }
  
  }
  
  if (type == "eisaR"){
    temp <- res_df[res_df$FDR < 0.05,]
    print(paste("Total genes included: ", nrow(temp_df), "   ",
                "FDR < 0.05: ", nrow(temp), "   ",
                "Up: ", sum(temp$logFC > 0), "   ",
                "Down: ", sum(temp$logFC < 0), "   ",
                "No change: ", sum(temp$logFC == 0), 
                sep = ""))
  }
}

# ----- load data -----
csvfile <- "sample_info_eisa.csv"
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = TRUE)

# R prefers first level of factor to be ref (ie. control, wild mode, untreated) sample. If colData table was assembled with treated samples set as reference, we need to use relevel like below. If not, skip this step
library("magrittr")
library("stringr")
coldata$cell %<>% relevel("Wt")
coldata$lps %<>% relevel("untrt")

coldata_Wt <- coldata[coldata$cell == "Wt", ]
coldata_H14 <- coldata[coldata$cell == "H14", ]
coldata_untrt <- coldata[coldata$lps == "untrt",]
coldata_trt <- coldata[coldata$lps == "trt", ]

# start here when changing sample_info
sample_info <- coldata_trt
filename <- "4 trt_H14 vs Wt"
list2env(load_count_data(sample_info, filename), .GlobalEnv)

# ----- EISA step by step -----
library("AnnotationDbi")
library("org.Mm.eg.db")
list2env(prepare_y(Rex, Rin, comp, ref, th=0), .GlobalEnv)

# 4. Statistical analysis in edgeR
library("edgeR")

#check if any discrepancies btw samples
cnt <- data.frame(Ex = Rex)
res_exon <- DE_analysis(cnt, quantGenes, filename, cnt_type = "exon", save_flag = TRUE)
print_summary(res_exon$res_edgeR, Dex, quantGenes, type="edgeR_exon")
plot_libsize(res_exon$y, cond, ref, comp, type="exon")

cnt <- data.frame(In = Rin)
res_intron <- DE_analysis(cnt, quantGenes, filename, cnt_type = "intron", save_flag = TRUE)
print_summary(res_intron$res_edgeR, Din, quantGenes, type="edgeR_intron")
plot_libsize(res_intron$y, cond, ref, comp, type="intron")

cnt <- data.frame(Ex = Rex, In = Rin)
res_EI <- DE_analysis(cnt, quantGenes, filename, cnt_type = "EI", save_flag = TRUE)

sigGenes <- res_exon$res_edgeR[res_exon$res_edgeR$FDR < 0.05,]
cnt <- data.frame(Ex = Rex[sigGenes[sigGenes$logFC > 0,]$ENTREZID,], In = Rin[sigGenes[sigGenes$logFC > 0,]$ENTREZID,])
res_EI_up <- DE_analysis(cnt, sigGenes[sigGenes$logFC > 0,]$ENTREZID, filename, cnt_type = "EI_up", save_flag = TRUE)
print_summary(res_EI_up$res_edgeR, Dex.Din, sigGenes[sigGenes$logFC > 0,]$ENTREZID, type="edgeR_EI_up")

cnt <- data.frame(Ex = Rex[sigGenes[sigGenes$logFC < 0,]$ENTREZID,], In = Rin[sigGenes[sigGenes$logFC < 0,]$ENTREZID,])
res_EI_down <- DE_analysis(cnt, sigGenes[sigGenes$logFC < 0,]$ENTREZID, filename, cnt_type = "EI_down", save_flag = TRUE)
print_summary(res_EI_down$res_edgeR, Dex.Din, sigGenes[sigGenes$logFC < 0,]$ENTREZID, type="edgeR_EI_down")

# ----- EISA one step -----
library("eisaR")
res <- runEISA(Rex[sigGenes$ENTREZID,], Rin[sigGenes$ENTREZID,], cond)
plotEISA(res)

res_eisaR <- as.data.frame(res$tab.ExIn[order(res$tab.ExIn$FDR),])
res_eisaR$symbol <- mapIds(org.Mm.eg.db, keys=rownames(res_eisaR),
                           column="SYMBOL", keytype="ENTREZID", multiVals="first")
res_eisaR$genename <- mapIds(org.Mm.eg.db, keys=rownames(res_eisaR),
                             column="GENENAME", keytype="ENTREZID", multiVals="first")
res_eisaR <- res_eisaR[, c(6, 7, 1, 2, 3, 4, 5)]
write.csv(res_eisaR, file=paste("results_eisa/", filename, "/results_eisaR.csv", sep=""), row.names=TRUE)

print_summary(res_eisaR, "", "", type="eisaR")


# 5. Visualise
# volcano plot
library("ggrepel")
options(ggrepel.max.overlaps = Inf)

# top n genes based on FDR to label on volcano plot
generate_mat <- function(res_edgeR, n = 20, th = 0.05){
  mat <- res_edgeR[,-c(1,3)]
  colnames(mat)[2] <- "L2FC"
  mat$topGenes <- mat$symbol
  mat$topGenes[-c(1:n)] <- ""
  sig <- mat$FDR < th
  print(paste(sum(sig), " significant genes", sep=""))
  mat$cols <- ifelse(sig, ifelse(mat$L2FC > 0, "up", "down"), paste("FDR > ", toString(th), sep=""))
  return(mat)
}

plot_volcano <- function(mat, type = 'exon', save_flag = FALSE, n=20){
  if (type == 'exon') {
    xlabel <- expression(paste("RNA change (log2 ",Delta,"exon)"))
  } else if (type == 'intron') {
    xlabel <- expression(paste("RNA change (log2 ",Delta,"intron)"))
  } else if (type == 'EI') {
    xlabel <- expression(paste("RNA change (log2 ",Delta,"exon/",Delta,"intron)"))
  }

  ggplot(mat, aes(x = L2FC, y = -log10(FDR), label = topGenes)) +
    geom_point(aes(colour=cols, size=logCPM)) +
    scale_size_continuous(range = c(0,4)) +
    scale_color_manual(values=alpha(c("#497AB3FF","#22222244","#E41A1CFF"), 0.4)) +
    ggtitle(str_sub(filename, 3,-1)) + 
    xlab(xlabel) + 
    ylab("Significance (-log10 FDR)") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel() 
  
  if (save_flag) {ggsave(paste("results_eisa/", filename, "/volcano_top", n, "_", type, ".png", sep=""), width=10, height=5)}
}

n=0
mat <- generate_mat(res_exon$res_edgeR, n=n)
plot_volcano(mat, type='exon', save_flag = TRUE, n=n)
mat <- generate_mat(res_intron$res_edgeR, n=n)
plot_volcano(mat, type='intron', save_flag = TRUE, n=n)
mat <- generate_mat(res_EI$res_edgeR, n=n)
plot_volcano(mat, type='EI', save_flag = TRUE, n=n)

n=20
mat <- generate_mat(res_exon$res_edgeR, n=n)
plot_volcano(mat, type='exon', save_flag = TRUE, n=n)
mat <- generate_mat(res_intron$res_edgeR, n=n)
plot_volcano(mat, type='intron', save_flag = TRUE, n=n)
mat <- generate_mat(res_EI$res_edgeR, n=n)
plot_volcano(mat, type='EI', save_flag = TRUE, n=n)


# Delta I vs. Delta E
generate_mat_DEDI <- function(res_col, res_label, n=20, th=0.05){
  label <- res_label[,c(1,2)]
  label$topGenes <- label$symbol
  label$topGenes[-c(1:n)] <- ""
  
  mat <- res_col[,-c(1,3)]
  colnames(mat)[2] <- "L2FC"
  sig <- mat$FDR < th
  print(paste(sum(sig), " significant genes", sep=""))
  mat$cols <- ifelse(sig, ifelse(mat$L2FC > 0, "up", "down"), paste("FDR > ", toString(th), sep="")) 
  
  mat_combined <- merge(mat, label[c("topGenes")], by="row.names", all=TRUE)
  rownames(mat_combined) <- mat_combined$Row.names
  mat_combined <- mat_combined[,-c(1)]
  return(mat_combined)
}

plot_DEDI <- function(eimat, label='exon', n=20){
  ggplot(eimat, aes(x = Din, y = Dex, label = topGenes)) +
    geom_point(aes(colour = cols, size = logCPM)) +
    scale_size_continuous(range = c(0,4)) +
    scale_color_manual(values= alpha(c("#497AB3FF", "#22222244", "#E41A1CFF"), 0.4)) +
    ggtitle(str_sub(filename, 3,-1)) + 
    xlab(expression(paste(Delta,"intron (log2 trt/untrt)", sep=""))) + 
    ylab(expression(paste(Delta,"exon (log2 trt/untrt)", sep=""))) + 
    geom_text_repel() 
  
  ggsave(paste("results_eisa/", filename, "/DIvDE_top", n, "_", label, ".png", sep=""), width=10, height=5)
}

n = 20
mat <- generate_mat_DEDI(res_EI$res_edgeR, res_exon$res_edgeR, n=n, th=0.05)
eimat <- data.frame(Din = rowMeans(Din), Dex = rowMeans(Dex), Dex.Din = rowMeans(Dex.Din))
eimat <- merge(eimat, mat[c("logCPM", "cols", "topGenes")], by="row.names", all=TRUE)
plot_DEDI(eimat, label='exon', n=n)

mat <- generate_mat_DEDI(res_EI$res_edgeR, res_intron$res_edgeR, n=n, th=0.05)
eimat <- data.frame(Din = rowMeans(Din), Dex = rowMeans(Dex), Dex.Din = rowMeans(Dex.Din))
eimat <- merge(eimat, mat[c("logCPM", "cols", "topGenes")], by="row.names", all=TRUE)
plot_DEDI(eimat, label='intron', n=n)

mat <- generate_mat_DEDI(res_EI$res_edgeR, res_EI$res_edgeR, n=n, th=0.05)
eimat <- data.frame(Din = rowMeans(Din), Dex = rowMeans(Dex), Dex.Din = rowMeans(Dex.Din))
eimat <- merge(eimat, mat[c("logCPM", "cols", "topGenes")], by="row.names", all=TRUE)
plot_DEDI(eimat, label='EI', n=n)

# no label
ggplot(eimat, aes(x = Din, y = Dex)) +
  geom_point(aes(colour = cols, size = logCPM)) +
  scale_size_continuous(range = c(0,4)) +
  scale_color_manual(values=alpha(c("#497AB3FF","#22222244","#E41A1CFF"),0.4)) +
  ggtitle(str_sub(filename, 3,-1)) + 
  xlab(expression(paste(Delta,"intron (log2 trt/untrt)", sep=""))) + 
  ylab(expression(paste(Delta,"exon (log2 trt/untrt)", sep="")))

ggsave(paste("results_eisa/", filename, "/DIvDE.png", sep=""), width=10, height=9)

sum(mat$cols=="up")
sum(mat$cols=="down")
cor.test(eimat$Dex, eimat$Din)

# only differentially expressed exons
mat <- generate_mat_exon(res_exon$res_edgeR, n=100)
eimat <- data.frame(Din = rowMeans(Din)[sigGenes$ENTREZID], Dex = rowMeans(Dex)[sigGenes$ENTREZID], Dex.Din = rowMeans(Dex.Din)[sigGenes$ENTREZID])
eimat <- merge(eimat, mat[c("logCPM", "cols", "topGenes")], by="row.names", all=TRUE)

# no colour
ggplot(eimat, aes(x = Din, y = Dex, label = topGenes)) +
  geom_point(aes(size = logCPM, color = alpha("black", 0.4))) +
  scale_size_continuous(range = c(0,4)) +
  ggtitle(str_sub(filename, 3,-1)) + 
  xlab(expression(paste(Delta,"intron (log2 trt/untrt)", sep=""))) + 
  ylab(expression(paste(Delta,"exon (log2 trt/untrt)", sep=""))) + 
  geom_text_repel() 


mat <- generate_mat(res_exon$res_edgeR, th=0.05)
eimat <- data.frame(Din = rowMeans(Din), Dex = rowMeans(Dex), Dex.Din = rowMeans(Dex.Din))
eimat <- merge(eimat, mat[c("logCPM", "cols", "topGenes")], by="row.names", all=TRUE)

#glimma XY plot - interactive HTML graphics
library(Glimma)
res_df <- res_exon$res_edgeR[res_exon$res_edgeR$FDR > 0,]
res_df[res_df$FDR == 0,]$FDR <- 1e-320

glXYPlot(res_df$logFC, -log10(res_df$FDR), xlab = "logFC", ylab = "-lg10(FDR)")

glXYPlot(x=fit$coef, y=fit$lod, xlab="logFC", ylab="logodds",
         status=dt, anno=anno, side.main="ID",
         counts=vm, groups=groups, sample.cols=sample.cols)

# check similarity between different packages
filepath <- "New folder"
res_quasR <- read.csv(paste(filepath, "/results_edgeR_exon.csv", sep=""), row.names=1, stringsAsFactors =  FALSE)
res_edgeR <- read.csv(paste(filepath, "/results_edgeR.csv", sep=""), row.names=1, stringsAsFactors =  FALSE)
res_combined <- read.csv(paste(filepath, "/results_combined.csv", sep=""), row.names=1, stringsAsFactors =  FALSE, skip=1)
res_deseq <- read.csv(paste(filepath, "/results_deseq2.csv", sep=""), row.names=1, stringsAsFactors =  FALSE)

gquasR <- res_quasR['symbol']
gsalmon <- res_salmon['symbol']
gcombined <- res_combined['symbol']
gdeseq <- res_deseq['symbol']

n <- 500
int <- Reduce(intersect, list(gquasR[1:n,1],gdeseq[1:n,1]))

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

# compare DEGs
filepaths = c('1 WT_trt vs untrt', '2 H14_trt vs untrt', '3 untrt_H14 vs WT', '4 trt_H14 vs WT')
res_1 <- read.csv(paste("results_eisa/", filepaths[2], "/results_edgeR_exon.csv", sep=""))
res_1 <- res_1[res_1$FDR < 0.05,]
res_2 <- read.csv(paste("results_eisa/", filepaths[1], "/results_edgeR_exon.csv", sep=""))
res_2 <- res_2[res_2$FDR < 0.05,]
res_3 <- read.csv(paste("results_eisa/", filepaths[3], "/results_edgeR_exon.csv", sep=""))
res_3 <- res_3[res_3$FDR < 0.05,]
res_4 <- read.csv(paste("results_eisa/", filepaths[4], "/results_edgeR_exon.csv", sep=""))
res_4 <- res_4[res_4$FDR < 0.05,]

compare_DEG <- function(res_1, res_2, name1 = "WT after Rx", name2 = "H14 after Rx", filename = "1 vs 2"){
  res1_up2 <- Reduce(intersect, list(res_1$symbol, res_2[res_2$logFC > 0,]$symbol))
  res1_down2 <- Reduce(intersect, list(res_1$symbol, res_2[res_2$logFC < 0,]$symbol))
  
  up1up2 <- Reduce(intersect, list(res_1[res_1$logFC > 0,]$symbol, res_2[res_2$logFC > 0,]$symbol))
  up1down2 <- Reduce(intersect, list(res_1[res_1$logFC > 0,]$symbol, res_2[res_2$logFC < 0,]$symbol))
  down1down2 <- Reduce(intersect, list(res_1[res_1$logFC < 0,]$symbol, res_2[res_2$logFC < 0,]$symbol))
  down1up2 <- Reduce(intersect, list(res_1[res_1$logFC < 0,]$symbol, res_2[res_2$logFC > 0,]$symbol))
  
  res1_notin2 <- setdiff(res_1$symbol, res_2$symbol)
  print(paste("up: ", length(res1_up2), "/", length(res_1$symbol[res_1$logFC > 0]),
              "  down: ", length(res1_down2), "/", length(res_1$symbol[res_1$logFC < 0]),
              "  not significant: ", length(res1_notin2), 
              sep=""))
  
  res_1$cols[res_1$symbol %in% res1_up2] <- 'up'
  res_1$cols[res_1$symbol %in% res1_down2] <- 'down'
  res_1$cols[res_1$symbol %in% res1_notin2] <- 'not significant'
  write.csv(res_1[res_1$symbol %in% up1up2,],file="results_eisa/up1up2.csv", row.names=TRUE)
  write.csv(res_1[res_1$symbol %in% up1down2,],file="results_eisa/up1down2.csv", row.names=TRUE)
  write.csv(res_1[res_1$symbol %in% down1up2,],file="results_eisa/down1up2.csv", row.names=TRUE)
  write.csv(res_1[res_1$symbol %in% down1down2,],file="results_eisa/down1down2.csv", row.names=TRUE)
  write.csv(res_1[res_1$symbol %in% res1_notin2,],file="results_eisa/1only.csv", row.names=TRUE)
  
  write.csv(res_1,file=paste("results_eisa/", filename, ".csv", sep=""), row.names=TRUE)
  
  ggplot(res_1, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(colour=cols, size=logCPM)) +
    scale_size_continuous(range = c(0,4)) +
    scale_color_manual(values=alpha(c("#497AB3FF","#22222244","#E41A1CFF"), 0.4)) +
    ggtitle(paste(name1, "_vs_", name2, sep="")) + 
    xlab(expression(paste("RNA change (log2 ",Delta,"exon)"))) + 
    ylab("Significance (-log10 FDR)") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")  

  ggsave(paste("results_eisa/volcano_", name1, "v", name2,".png", sep=""), width=10, height=5)
}

compare_DEG(res_1, res_2, name1 = "WT after Rx", name2 = "H14 after Rx", filename = "1 vs 2")
compare_DEG(res_2, res_1, name1 = "H14 after Rx", name2 = "WT after Rx", filename = "2 vs 1")
