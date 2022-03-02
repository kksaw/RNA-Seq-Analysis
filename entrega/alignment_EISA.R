suppressPackageStartupMessages({
  library(argparser)
  library(glue)
  library(dplyr)
  library(eisaR)
  library(QuasR)

  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# ----- parser -----
p <- arg_parser("QuasR alignment and counting for EISA")
# p <- add_argument(p, "-w", help="work directory")
p <- add_argument(p, "-t", help="no. of threads")
# p <- add_arguemnt(p, "-f", help="path to .fa genome") #?
# p <- add_argument(p, "-g", help="path to .gtf file") #txdb
p <- add_argument(p, "-q", help="path to QuasR sample file")
args <- parse_args(p)

# setwd(args$w)

# ----- load data -----
# gse <- tail(strsplit(getwd(), "/")[[1]], n=1)
genomeFile <- "BSgenome.Mmusculus.UCSC.mm10"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# genomeFile <- args$f
# txdb <- makeTxDbFromGFF(args$g)
sampleFile <- args$q

# make cluster
cl <- makeCluster(as.numeric(args$t))

# ----- alignment -----
proj <- qAlign(sampleFile = sampleFile,
               genome = genomeFile,
               aligner = "Rhisat2", splicedAlignment = TRUE,
               alignmentsDir = "./bam", cacheDir = "./cache", clObj = cl)

# save alignment stats
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE)

# ----- count exons and introns
regions = getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)

# count alignments in exons and gene bodies
cntEx <- qCount(proj, regions$exons,
                orientation = "same", useRead = "first", clObj = cl)
cntGb <- qCount(proj, regions$genebodies,
                orientation = "same", useRead = "first", clObj = cl)
cntIn <- cntGb - cntEx  # intron counts = difference btw gene bodies and exons

Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# save to output
write.table(Rex, file=glue("quants_eisa/ExonicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(Rin, file=glue("quants_eisa/IntronicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# delete cache
system("rm -r cache/*")
