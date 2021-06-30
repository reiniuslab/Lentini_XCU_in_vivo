library(data.table)

##

loom2matrix <- function(x){
  require(loomR)
  lfile <- connect(x)
  out <- t(lfile$matrix[,])
  colnames(out) <- lfile$col.attrs$CellID[]
  row.names(out) <- lfile$row.attrs$Gene[]
  lfile$close_all()
  return(out)
}

cbind.intersect <- function(x,y){
  idx <- intersect(row.names(x), row.names(y))
  cbind(x[idx,], y[idx,])
}

get.x.frac.ref <- function(ref, alt, chrom, cutoff=0, prob=0.9){
  # get chrX fraction after excluding the top {prob} expressed genes
  chrx <- which(chrom == "chrX")
  total <- ref+alt
  x.avg <- rowMeans(total[chrx,], na.rm = T)
  idx <- x.avg < quantile(x.avg[x.avg > cutoff],probs = prob)
  frac <- colSums(ref[chrx,][idx,])/colSums(total[chrx,][idx,])
  frac
}

get.x.status <- function(frac){
  x <- NA
  x[frac >= 0.9] <- "XaXi"
  x[frac <= 0.1] <- "XiXa"
  x[frac < 0.9 & frac >= 0.6] <- "XaXs"
  x[frac > 0.1 & frac <= 0.4] <- "XsXa"
  x[frac > 0.4 & frac < 0.6] <- "XaXa"
  x
}

fpkm2tpm <- function(fpkm){
  (fpkm / sum(fpkm,na.rm = T)) * 1e6
}

getAllelicRatio <- function(ArchRProj, useMatrix, ... ){
  se <- getMatrixFromProject(ArchRProj, useMatrix, logFile = NULL, ... )
  se <- se[,order(colnames(se))]
  nlen <- ncol(se)
  refratio <- assay(se)[,1:(nlen/2)] / (assay(se)[,1:(nlen/2)] + assay(se)[,((nlen/2)+1):nlen])
  
  out <- SummarizedExperiment(assays = list(RefRatio = refratio), rowData = rowData(se), colData = colData(se)[1:(nlen/2),])
  colnames(out) <- gsub(".*#", "", colnames(out))
  return(out)
}

getGroupedMatrix <- function(x, grp, fun, ... ){
  # a lot faster than getGroupSE() and apply(x,1,tapply) but requires row* / col* functions from the Matrix package
  grp <- factor(grp)
  sapply(levels(grp), function(i) fun(x[, which(grp == i)], ...) )
}

getGroupedMatrix_OLD <- function(x, grp, fun, ... ){
  # slower but works with any function
  t(apply(x, 1, function(y) tapply(y, grp, fun, ...) ))
}

sce2se <- function(sce, exprs_values = "counts", features, addchr = T){
  library(SummarizedExperiment)
  library(GenomicRanges)
  rr <- GRanges(features)
  
  if(addchr){
    rr <- renameSeqlevels(rr, paste0("chr", seqlevels(rr)))
  }
  se <- SummarizedExperiment(assays = SimpleList(counts = Matrix(round(assay(sce, exprs_values)), sparse = T)), rowRanges = rr, colData = colData(sce))
  rownames(se) <- features$gene_name
  
  return(se)
}

##

## scATAC
# load data
library(ArchR)
ap.total.filt <- loadArchRProject("scATAC/XupregTotal")
ap.allele.filt <- loadArchRProject("scATAC/XupregAllele")

## scRNA
# load data
library(rtracklayer)
gtf <- import("smartseq3/Mus_musculus.GRCm38.97.chr.gtf.gz")
gene.anno.ss3 <- unique(with(gtf, data.table(gene_id, gene_name, chrom = as.character(seqnames), start, end)))
gene.anno.ss3 %<>% .[chrom %in% c(1:19,"X","Y")]

rm(list = c("gtf"))

# load meta
meta.paired.1 <- fread("scATAC/paired_SS3_data/scRNA_metadata_d1d2_fixed.tsv")
meta.paired.1[, Barcode := paste0(Index1, Index2)]
meta.paired.1[, ATAC := gsub("R_", "D_", Sample_ID, fixed = T) %in% ap.total.filt$Names]

meta.paired.2 <- fread("scATAC/paired_SS3_data/scRNA_metadata_d4d7_fixed.tsv")
meta.paired.2[, Barcode := paste0(Index1, Index2)]
meta.paired.2[, ATAC := gsub("R_", "D_", Sample_ID, fixed = T) %in% ap.total.filt$Names]

# load total counts as SE
library(scater)
library(scran)
library(magrittr)
sce.paired.1 <- SingleCellExperiment(
  assays = list(
    counts = loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d1d2.readcount.exon.all.loom"),
    umi = loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d1d2.umicount.exon.all.loom")
  ),
  colData = DataFrame(meta.paired.1[match(colnames(loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d1d2.readcount.exon.all.loom")), Barcode)])
)
colnames(sce.paired.1) <- sce.paired.1$Sample_ID

sce.paired.2 <- SingleCellExperiment(
  assays = list(
    counts = loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d4d7.readcount.exon.all.loom"),
    umi = loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d4d7.umicount.exon.all.loom")
  ),
  colData = DataFrame(meta.paired.2[match(colnames(loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d4d7.readcount.exon.all.loom")), Barcode)])
)
colnames(sce.paired.2) <- sce.paired.2$Sample_ID

# merge
sce.paired <- cbind(sce.paired.1[intersect(rownames(sce.paired.1), rownames(sce.paired.2))], sce.paired.2[intersect(rownames(sce.paired.1), rownames(sce.paired.2))])
rowData(sce.paired) <- gene.anno.ss3[match(rownames(sce.paired), gene_id), list(gene_id, gene_name, chrom = paste0("chr", chrom))]
sce.paired %<>% .[!rowSums(is.na(rowData(sce.paired))) > 0]
sce.paired %<>% calculateQCMetrics()
sce.paired$exclude <- isOutlier(sce.paired$total_counts, nmads = 3, type = "lower", log=T, batch = sce.paired$Plate) | isOutlier(sce.paired$total_features_by_counts, nmads = 3, type = "lower", log=T, batch = sce.paired$Plate)

sce.filt.paired <- subset(sce.paired,, !exclude)

meta.paired <- rbind(meta.paired.1, meta.paired.2)
meta.paired %<>% .[match(colnames(sce.paired), Sample_ID)]

rm(list=c("meta.paired.1", "meta.paired.2", "sce.paired.1", "sce.paired.2"))

# load allelic counts
counts.paired.1 <- list("c57" = read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d1d2.BL6_reads.txt",row.names=1), "cast" = read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d1d2.CAST_reads.txt",row.names=1))
colnames(counts.paired.1[[1]]) <- meta.paired[Batch == "A"][match(colnames(counts.paired.1[[1]]), Barcode), Sample_ID]; colnames(counts.paired.1[[2]]) <- meta.paired[Batch == "A"][match(colnames(counts.paired.1[[2]]), Barcode), Sample_ID]
counts.paired.2 <- list("c57" = read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d4d7.BL6_reads.txt",row.names=1), "cast" = read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d4d7.CAST_reads.txt",row.names=1))
colnames(counts.paired.2[[1]]) <- meta.paired[Batch == "B"][match(colnames(counts.paired.2[[1]]), Barcode), Sample_ID]; colnames(counts.paired.2[[2]]) <- meta.paired[Batch == "B"][match(colnames(counts.paired.2[[2]]), Barcode), Sample_ID]

counts.paired <- list(
  "c57" = cbind.intersect(counts.paired.1$c57, counts.paired.2$c57),
  "cast" = cbind.intersect(counts.paired.1$cast, counts.paired.2$cast)
)

rm(list = c("counts.paired.1", "counts.paired.2"))

# load rpkm
rpkm.all.paired.1 <- loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d1d2.rpkm.exon.all.loom")
colnames(rpkm.all.paired.1) <- meta.paired[Batch == "A"][match(colnames(rpkm.all.paired.1), Barcode), Sample_ID]
rpkm.all.paired.2 <- loom2matrix("scATAC/paired_SS3_data/zUMIs_output/expression/scSS3_of_ATAC_plates_d4d7.rpkm.exon.all.loom")
colnames(rpkm.all.paired.2) <- meta.paired[Batch == "B"][match(colnames(rpkm.all.paired.2), Barcode), Sample_ID]

tpm.all.paired <- apply(cbind.intersect(rpkm.all.paired.1, rpkm.all.paired.2),2,fpkm2tpm)

rm(list = c("rpkm.all.paired.1", "rpkm.all.paired.2"))

# load gene-wise allele ratios and remove genes only mapping to one allele
altratio.paired.1 <- read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d1d2.fract_CAST_reads.txt", row.names=1)
colnames(altratio.paired.1) <- meta.paired[Batch == "A"][match(colnames(altratio.paired.1), Barcode), Sample_ID]
altratio.paired.2 <- read.delim("scATAC/paired_SS3_data/zUMIs_output/allelic/scSS3_of_ATAC_plates_d4d7.fract_CAST_reads.txt", row.names=1)
colnames(altratio.paired.2) <- meta.paired[Batch == "B"][match(colnames(altratio.paired.2), Barcode), Sample_ID]

altratio.paired <- cbind.intersect(altratio.paired.1, altratio.paired.2)
altratio.paired.filt <- altratio.paired[which(apply(altratio.paired,1,var, na.rm=T)>0),]

rm(list = c("altratio.paired.1", "altratio.paired.2", "altratio.paired"))

# scale data by allele counts
tpm.allele.paired <- list(
  "c57" = as.matrix(tpm.all.paired[intersect(row.names(altratio.paired.filt), row.names(tpm.all.paired)), colnames(altratio.paired.filt)] * (1 - altratio.paired.filt)),
  "cast" = as.matrix(tpm.all.paired[intersect(row.names(altratio.paired.filt), row.names(tpm.all.paired)), colnames(altratio.paired.filt)] * altratio.paired.filt)
)

# add X ratio
c57.x.frac.paired <- get.x.frac.ref(counts.paired$c57, counts.paired$cast, gene.anno.ss3[match(row.names(counts.paired$c57), gene_id), paste0("chr",chrom)], prob = 0.99)
meta.paired[, c57.x.frac := c57.x.frac.paired[match(meta.paired[,Sample_ID], names(c57.x.frac.paired))]]
meta.paired[, x.status := get.x.status(c57.x.frac)]
meta.paired[!is.na(x.status), Group := paste(Sex, x.status, sep = "_")]
meta.paired[Group %in% names(table(Group))[table(Group) < 10], Group := NA] # Exclude too small groups
meta.paired[as.data.table(colData(sce.paired)), exclude := exclude, on = "Sample_ID == Sample_ID"]

rm(list = c("c57.x.frac.paired"))

# melt data
tpm.melt.paired <- na.omit(as.data.table(melt(tpm.allele.paired)))
colnames(tpm.melt.paired) <- c("gene_id", "sample", "value", "allele")
tpm.melt.paired %<>% .[meta.paired, on = c("sample == Sample_ID")]
tpm.melt.paired[gene.anno.ss3, c("gene_name", "chrom") := list(gene_name, chrom), on = c("gene_id == gene_id")]
colnames(tpm.melt.paired) %<>% tolower()
tpm.melt.paired[,expressed := mean(value,na.rm = T)>0, by = c("gene_id","diff_day")]

# Update ArchR objects
if(length(dir("scATAC/XupregAllele/GroupCoverages/Group/")) == 0){
  ## add multi-modal information
  se.paired <- sce2se(sce.filt.paired, "umi", gene.anno.ss3[match(rownames(sce.filt.paired), gene_id)])
  se.paired %<>% .[,gsub("R_", "_", colnames(se.paired)) %in% gsub("D_", "_", ap.total.filt$Names)]
  colnames(se.paired) <- ap.total.filt$cellNames[match(gsub("R_", "_", colnames(se.paired)), gsub("D_", "_", ap.total.filt$Names))]
  
  ap.total.filt2 <- addGeneExpressionMatrix(ap.total.filt, seRNA = se.paired, threads = 1)
  ap.total.filt2 %<>% .[!is.na(ap.total.filt2$Gex_nUMI)] # only cells that pass paired qc
  
  # add dimensionality reductions
  ap.total.filt2 %<>% addIterativeLSI(useMatrix = "TileMatrix", saveIterations = FALSE, depthCol = "nFrags", excludeChr = c("chrY", "chrMT"), name = "LSI_ATAC", force = T)
  ap.total.filt2 %<>% addIterativeLSI(useMatrix = "GeneExpressionMatrix", saveIterations = FALSE, depthCol = "Gex_nUMI", excludeChr = c("chrY", "chrMT"), varFeatures = 1000, firstSelection = "Var", selectionMethod = "vmr", binarize = FALSE, name = "LSI_RNA", force = T)
  ap.total.filt2 %<>% addCombinedDims(reducedDims = c("LSI_ATAC", "LSI_RNA"), scaleDims = F, dimsToUse = 1:5, corCutOff = 0.4, name = "LSI_Combined")
  
  ## add X ratio to scATAC
  ap.total.filt$c57.x.frac <- meta.paired[match(ap.total.filt$Names, gsub("R","D",Sample_ID)), c57.x.frac]
  ap.total.filt$x.status <- get.x.status(ap.total.filt$c57.x.frac)
  ap.total.filt$Group <- NA
  ap.total.filt$Group[!is.na(ap.total.filt$x.status)] <- with(getCellColData(ap.total.filt[!is.na(ap.total.filt$x.status)]), paste(Sex, x.status, sep="_"))
  ap.total.filt$Group[with(getCellColData(ap.total.filt), Group %in% names(table(Group))[table(Group)<10])] <- NA # Exclude too small groups

  ap.total.filt2$c57.x.frac <- meta.paired[match(ap.total.filt2$Names, gsub("R","D",Sample_ID)), c57.x.frac]
  ap.total.filt2$x.status <- get.x.status(ap.total.filt2$c57.x.frac)
  ap.total.filt2$Group <- NA
  ap.total.filt2$Group[!is.na(ap.total.filt2$x.status)] <- with(getCellColData(ap.total.filt2[!is.na(ap.total.filt2$x.status)]), paste(Sex, x.status, sep="_"))
  ap.total.filt2$Group[with(getCellColData(ap.total.filt2), Group %in% names(table(Group))[table(Group)<10])] <- NA # Exclude too small groups
  
  ap.allele.filt$c57.x.frac <- meta.paired[match(ap.allele.filt$Names, gsub("R","D",Sample_ID)), c57.x.frac]
  ap.allele.filt$x.status <- get.x.status(ap.allele.filt$c57.x.frac)
  ap.allele.filt$Group <- NA
  ap.allele.filt$Group[!is.na(ap.allele.filt$x.status)] <- with(getCellColData(ap.allele.filt[!is.na(ap.allele.filt$x.status)]), paste(Sex, x.status, Allele, sep="_"))
  ap.allele.filt$Group[with(getCellColData(ap.allele.filt), Group %in% names(table(Group))[table(Group)<10])] <- NA # Exclude too small groups
  
  # add dim reduction
  ap.total.filt2 %<>% addUMAP(reducedDims = "LSI_ATAC", name = "UMAP_ATAC", nNeighbors = 15, threads = 1, force=T)
  ap.total.filt2 %<>% addUMAP(reducedDims = "LSI_RNA", name = "UMAP_RNA", nNeighbors = 15, threads = 1, force=T)
  ap.total.filt2 %<>% addUMAP(reducedDims = "LSI_Combined", name = "UMAP_Combined", nNeighbors = 15, threads = 1, force = T)
  
  ## add group coverage
  ap.total.filt %<>% addGroupCoverages(groupBy = "Group")
  ap.total.filt2 %<>% addGroupCoverages(groupBy = "Group")
  ap.allele.filt %<>% addGroupCoverages(groupBy = "Group")
  
  ## Save updated projects
  saveArchRProject(ArchRProj = ap.total.filt, outputDirectory = "scATAC/XupregTotal", load = FALSE)
  saveArchRProject(ArchRProj = ap.total.filt2, outputDirectory = "scATAC/XupregTotalIntegrated", load = FALSE)
  saveArchRProject(ArchRProj = ap.allele.filt, outputDirectory = "scATAC/XupregAllele", load = FALSE)
}

ap.total.filt2 <- loadArchRProject("scATAC/XupregTotalIntegrated", showLogo = F)