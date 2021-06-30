library(ArchR)
library(magrittr)

addArchRGenome("mm10")

## Read bam files
bfls <- c("bam/merge/scATAC_all.sorted.bam")
#bfls <- c("scATAC/bam/merge/scATAC_all.sorted.bam")
names(bfls) <- "scATAC_Total"

bfls.allelic <- c("bam/merge/scATAC_all.genome1.sorted.bam","bam/merge/scATAC_all.genome2.sorted.bam")
#bfls.allelic <- c("scATAC/bam/merge/scATAC_all.genome1.sorted.bam","scATAC/bam/merge/scATAC_all.genome2.sorted.bam")
names(bfls.allelic) <- c("scATAC_C57","scATAC_CAST")

bamFlag <- list(isPaired = T, isProperPair = T,
                isDuplicate = F, isSecondaryAlignment = F,
                isNotPassingQualityControls = F, isUnmappedQuery = F, 
                hasUnmappedMate = F)

ArrowFiles <- createArrowFiles(
  inputFiles = bfls,
  gsubExpression = ":.*",
  excludeChr = c("chrMT", "chrY"),
  bamFlag = bamFlag,
  verbose = F
)

ArrowFiles.allelic <- createArrowFiles(
  inputFiles = bfls.allelic,
  gsubExpression = ":.*",
  minTSS = 0,
  minFrags = 0,
  excludeChr = c("chrMT", "chrY"),
  bamFlag = bamFlag,
  addGeneScoreMat = F,
  addTileMat = F,
  verbose = F
)

## Calculate doublets
addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

## Create project
ap.total <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "XupregTotal",  copyArrows = TRUE)
#ap.total <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "scATAC/XupregTotal",  copyArrows = TRUE)
ap.allele <- ArchRProject(ArrowFiles = ArrowFiles.allelic, outputDirectory = "XupregAllele",  copyArrows = TRUE)
#ap.allele <- ArchRProject(ArrowFiles = ArrowFiles.allelic, outputDirectory = "scATAC/XupregAllele",  copyArrows = TRUE)

## Add metadata
# Read metadata
library(data.table)
meta <- fread("scATAC_metadata_fixed.tsv")
#meta <- fread("scATAC/scATAC_metadata_fixed.tsv")

# total
ap.total$Names <- gsub(".*#|_S[[:digit:]]{1,4}","",ap.total$cellNames)
ap.total$Day <- meta[match(ap.total$Names, Sample_ID),paste0("D",Diff_day)]
ap.total$Sex <- meta[match(ap.total$Names, Sample_ID),Sex]

## Filter data
ap.total.filt <- ap.total[ap.total$PassQC == T]
ap.total.filt %<>% filterDoublets()

# allelic
ap.allele$Names <- gsub(".*#|_S[[:digit:]]{1,4}","",ap.allele$cellNames)
ap.allele$Allele <- gsub("scATAC_|#.*","",ap.allele$cellNames)
ap.allele$Day <- meta[match(ap.allele$Names, Sample_ID),paste0("D",Diff_day)]
ap.allele$Sex <- meta[match(ap.allele$Names, Sample_ID),Sex]

ap.allele.filt <- ap.allele[ap.allele$Names %in% ap.total.filt$Names]

# Add count matrices
ap.allele.filt %<>% addTileMatrix(force=T)
ap.allele.filt %<>% addGeneScoreMatrix()
ap.allele.filt %<>% addFeatureMatrix(features = keepSeqlevels(getGenes(ap.allele.filt), paste0("chr",c(1:19,"X")), pruning.mode = "coarse")[-11849], matrixName="GeneCountMatrix") # exclude chrY, chrM and 1 NA gene

## Call peaks
ap.total.filt %<>% addGroupCoverages(groupBy = "Day")
ap.total.filt %<>% addReproduciblePeakSet(groupBy = "Day", peakMethod = "Macs2")
ap.total.filt %<>% addPeakMatrix()

ap.allele.filt %<>% addPeakSet(peakSet=getPeakSet(ap.total.filt))
ap.allele.filt %<>% addPeakMatrix()

## Save output projects
saveArchRProject(ArchRProj = ap.total.filt, outputDirectory = "XupregTotal", load = FALSE)
#saveArchRProject(ArchRProj = ap.total.filt, outputDirectory = "scATAC/XupregTotal", load = FALSE)
saveArchRProject(ArchRProj = ap.allele.filt, outputDirectory = "XupregAllele", load = FALSE)
#saveArchRProject(ArchRProj = ap.allele.filt, outputDirectory = "scATAC/XupregAllele", load = FALSE)