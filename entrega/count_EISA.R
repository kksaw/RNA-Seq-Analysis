# ----- make qProject class from bam files -----
suppressPackageStartupMessages({
  library(argparser)
  library(glue)
  library(dplyr)
  library(eisaR)
  library(QuasR)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# ----- parser -----
p <- arg_parser("QuasR counting for EISA")
# p <- add_argument(p, "-w", help="work directory")
p <- add_argument(p, "-t", help="no. of threads")
# p <- add_arguemnt(p, "-f", help="path to .fa genome") #?
# p <- add_argument(p, "-g", help="path to .gtf file") #txdb
p <- add_argument(p, "-d", help="path to QuasR project file")
p <- add_argument(p, "-q", help="path to QuasR sample file")
args <- parse_args(p)

cl <- makeCluster(as.numeric(args$t))

qProject <- setClass("qProject", slots=list(
  reads = "data.frame",
  reads_md5subsum = "data.frame",
  alignments = "data.frame",
  samplesFormat = "character",
  genome = "character",
  genomeFormat = "character",
  aux = "data.frame",
  auxAlignments = "data.frame",
  aligner = "character",
  maxHits = "numeric",
  paired = "character",
  splicedAlignment = "logical" ,
  snpFile = "character",
  bisulfite = "character",
  alignmentParameter = "character",
  projectName = "character",
  alignmentsDir = "character",
  lib.loc = "character",
  cacheDir = "character",
  alnModeID = "character",
  geneAnnotation = "character",
  geneAnnotationFormat = "character"))

projdetails <- read.csv(args$d, colClasses=c("phred"="character"))
sampledetails <- read.delim(args$q)

proj <- qProject(
  reads = cbind(sampledetails, projdetails[3]),
  reads_md5subsum = projdetails[4],
  alignments = projdetails[1:2],
  samplesFormat = "fastq",
  genome = "BSgenome.Mmusculus.UCSC.mm10",
  genomeFormat = "BSgenome",
  aligner = "Rhisat2",
  maxHits = 1,
  paired = "no",
  splicedAlignment = TRUE,
  snpFile = NA_character_,
  bisulfite = "no",
  alignmentParameter = "-k 2",
  projectName = "qProject",
  alignmentsDir = "./bam",
  lib.loc = NA_character_,
  cacheDir = "./cache",
  alnModeID = "Rhisat2",
  geneAnnotation = NA_character_,
  geneAnnotationFormat = NA_character_)

genomeFile <- "BSgenome.Mmusculus.UCSC.mm10"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
regions = getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)

# count alignments in exons and gene bodies
# SAME ORIENTATION, FIRST READ
cntEx <- qCount(proj, regions$exons,
                orientation = "same", useRead = "first", clObj=cl)
cntGb <- qCount(proj, regions$genebodies,
                orientation = "same", useRead = "first", clObj=cl)
cntIn <- cntGb - cntEx  # intron counts = difference btw gene bodies and exons

Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# save to output
write.table(Rex, file=glue("quants_eisa/ExonicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(Rin, file=glue("quants_eisa/IntronicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)


# ANY ORIENTATION, FIRST READ
cntEx <- qCount(proj, regions$exons,
                orientation = "any", useRead = "first", clObj=cl)
cntGb <- qCount(proj, regions$genebodies,
                orientation = "any", useRead = "first", clObj=cl)
cntIn <- cntGb - cntEx  # intron counts = difference btw gene bodies and exons

Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# save to output
write.table(Rex, file=glue("quants_eisa/ExonicCounts_any.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(Rin, file=glue("quants_eisa/IntronicCounts_any.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)


# DEFAULT: ANY ORIENTATION, ANY READ
cntEx <- qCount(proj, regions$exons,
                orientation = "any", useRead = "any", clObj=cl)
cntGb <- qCount(proj, regions$genebodies,
                orientation = "any", useRead = "any", clObj=cl)
cntIn <- cntGb - cntEx  # intron counts = difference btw gene bodies and exons

Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# save to output
write.table(Rex, file=glue("quants_eisa/ExonicCounts_default.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(Rin, file=glue("quants_eisa/IntronicCounts_default.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

