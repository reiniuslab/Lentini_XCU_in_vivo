library(data.table)
library(magrittr)
library(biomaRt)

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

fpkm2tpm <- function(fpkm){
  (fpkm / sum(fpkm,na.rm = T)) * 1e6
}

umi2rel <- function(umi){
  (umi / sum(umi,na.rm = T)) * 1e4
}

add.info.ss3 <- function(x, coldata, rowdata = gene.anno.ss3){
  library(data.table)
  x[rowdata, c("gene","chr") := list(gene_name, paste0("chr",chrom)), on="Var1 == gene_id"]
  x[chr != "chrY", chrx := chr == "chrX"]
  x[coldata, c("sample_id","exclude","clone","sex","day","c57.x.frac","x.status") := list(sample_id,exclude,clone,sex,day,c57.x.frac,x.status),on = "Var2 == sample_bc"]
  names(x) %<>% tolower()
  return(x)
}

get.x.frac.ref <- function(ref, alt, chrom, cutoff=0, prob=0.9){
  # get chrX fraction after excluding the top {prob} expressed genes
  chrx <- which(chrom == "chrX")
  total <- ref+alt
  x.avg <- rowMeans(total[chrx,])
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

scale.by.allele <- function(x, ref, alt, labs=c("c57","cast")){
  total <- ref+alt
  refratio <- ref/total
  x.ref <- x*refratio
  x.alt <- x*(1-refratio)
  x.allele <- list(x.ref, x.alt)
  names(x.allele) <- labs
  return(x.allele)
}

load.umi.ss3 <- function(x,y, ... ){
  require(data.table)
  ls <- list(
    "c57" = as.matrix(read.delim(x,row.names = 1)),
    "cast" = as.matrix(read.delim(y,row.names = 1))
  )
  dat.melt <- as.data.table(melt(ls))
  dat.melt %<>% add.info.ss3( ... )
  dat.melt[, umi_rel := umi2rel(value), by="var2"] # median(colSums(Reduce("+",ls),na.rm = T))
  return(dat.melt)
}

median_cl_boot <- function(x, conf = 0.95) {
  # function from http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}

##

# load annotations
library(rtracklayer)
gtf <- import("smartseq3/Mus_musculus.GRCm38.97.chr.gtf.gz")
gene.anno.ss3 <- unique(with(gtf, data.table(gene_id, gene_name, chrom = as.character(seqnames), start)))
gene.anno.ss3 %<>% .[chrom %in% c(1:19,"X","Y")]

# load counts
counts.ss3.day0 <- list(
  "c57" = read.delim("smartseq3/results/zUMIs_output/allelic/mESC_day0.BL6_reads.txt",row.names=1),
  "cast" = read.delim("smartseq3/results/zUMIs_output/allelic/mESC_day0.CAST_reads.txt",row.names=1)
)

counts.ss3.cl5diff <- list(
  "c57" = read.delim("smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.BL6_reads.txt",row.names=1),
  "cast" = read.delim("smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.CAST_reads.txt",row.names=1)
)

# update metadata with X allele fractions
meta.ss3.day0 <- fread("smartseq3/metadata_ss3_day0.tsv")
meta.ss3.day0[as.data.table(get.x.frac.ref(counts.ss3.day0$c57, counts.ss3.day0$cast, gene.anno.ss3[match(row.names(counts.ss3.day0$c57), gene_id), paste0("chr",chrom)], prob = 0.99),keep.rownames = T), c57.x.frac := V2, on="sample_bc == V1"] # this data is more biased, use higher cut-off to correctly assign males
meta.ss3.day0[,x.status := get.x.status(c57.x.frac)]
meta.ss3.day0[fread("smartseq3/samples_exclude.tsv"),exclude:=exclude,on = "sample_id==sample_id"]

meta.ss3.cl5diff <- fread("smartseq3/metadata_ss3_clone5diff.tsv")
meta.ss3.cl5diff[as.data.table(get.x.frac.ref(counts.ss3.cl5diff$c57, counts.ss3.cl5diff$cast, gene.anno.ss3[match(row.names(counts.ss3.cl5diff$c57), gene_id), paste0("chr",chrom)]),keep.rownames = T), c57.x.frac := V2, on="sample_bc == V1"]
meta.ss3.cl5diff[,x.status := get.x.status(c57.x.frac)]
meta.ss3.cl5diff[fread("smartseq3/samples_exclude.tsv"),exclude:=exclude,on = "sample_id==sample_id"]

meta.ss3 <- rbindlist(list(meta.ss3.day0,meta.ss3.cl5diff))
# fwrite(meta.ss3, "smartseq3/metadata_ss3_full.tsv", quote = F, sep = "\t") # export updated metadata

# load RPKM data
rpkm.all.ss3.day0 <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_day0.rpkm.exon.all.loom")
tpm.all.ss3.day0 <- apply(rpkm.all.ss3.day0,2,fpkm2tpm)
tpm.melt.all.ss3.day0 <- as.data.table(melt(tpm.all.ss3.day0))
tpm.melt.all.ss3.day0 %<>% add.info.ss3(coldata = meta.ss3.day0)
tpm.melt.all.ss3.day0 %<>% .[exclude == F] # exclude low read depth cells

rpkm.all.ss3.cl5diff <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_diff_clone5.rpkm.exon.all.loom")
tpm.all.ss3.cl5diff <- apply(rpkm.all.ss3.cl5diff,2,fpkm2tpm)
tpm.melt.all.ss3.cl5diff <- as.data.table(melt(tpm.all.ss3.cl5diff))
tpm.melt.all.ss3.cl5diff %<>% add.info.ss3(coldata = meta.ss3.cl5diff)
tpm.melt.all.ss3.cl5diff %<>% .[exclude == F] # exclude low read depth cells

tpm.melt.all.ss3 <- rbindlist(list(
  tpm.melt.all.ss3.day0,
  tpm.melt.all.ss3.cl5diff
))

tpm.melt.all.ss3[,expressed := mean(value,na.rm = T)>0, by = c("gene","day")]

# load gene-wise allele ratios and remove genes only mapping to one allele
altratio.ss3.day0 <- read.delim("smartseq3/results/zUMIs_output/allelic/mESC_day0.fract_CAST_reads.txt", row.names=1)
altratio.ss3.day0.filt <- altratio.ss3.day0[which(apply(altratio.ss3.day0,1,var, na.rm=T)>0),]
# scale data by allele counts
tpm.allele.ss3.day0 <- list(
  "c57" = as.matrix(tpm.all.ss3.day0[intersect(row.names(altratio.ss3.day0.filt), row.names(tpm.all.ss3.day0)), colnames(altratio.ss3.day0.filt)] * (1 - altratio.ss3.day0.filt)),
  "cast" = as.matrix(tpm.all.ss3.day0[intersect(row.names(altratio.ss3.day0.filt), row.names(tpm.all.ss3.day0)), colnames(altratio.ss3.day0.filt)] * altratio.ss3.day0.filt)
)
tpm.melt.ss3.day0 <- as.data.table(melt(tpm.allele.ss3.day0))
tpm.melt.ss3.day0 %<>% add.info.ss3(coldata=meta.ss3.day0)
tpm.melt.ss3.day0 %<>% .[exclude == F]

altratio.ss3.cl5diff <- read.delim("smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.fract_CAST_reads.txt", row.names=1)
altratio.ss3.cl5diff.filt <- altratio.ss3.cl5diff[which(apply(altratio.ss3.cl5diff,1,var, na.rm=T)>0),]
tpm.allele.ss3.cl5diff <- list(
  "c57" = as.matrix(tpm.all.ss3.cl5diff[intersect(row.names(altratio.ss3.cl5diff.filt), row.names(tpm.all.ss3.cl5diff)), colnames(altratio.ss3.cl5diff.filt)] * (1 - altratio.ss3.cl5diff.filt)),
  "cast" = as.matrix(tpm.all.ss3.cl5diff[intersect(row.names(altratio.ss3.cl5diff.filt), row.names(tpm.all.ss3.cl5diff)), colnames(altratio.ss3.cl5diff.filt)] * altratio.ss3.cl5diff.filt)
)
tpm.melt.ss3.cl5diff <- as.data.table(melt(tpm.allele.ss3.cl5diff))
tpm.melt.ss3.cl5diff %<>% add.info.ss3(coldata=meta.ss3.cl5diff)
tpm.melt.ss3.cl5diff %<>% .[exclude == F]

tpm.melt.ss3 <- rbindlist(list(tpm.melt.ss3.day0,tpm.melt.ss3.cl5diff))
tpm.melt.ss3[,expressed := tpm.melt.all.ss3$expressed[match(paste0(gene,day),with(tpm.melt.all.ss3,paste0(gene,day)))]]

# load UMI data
umi.all.ss3.day0 <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_day0.umicount.exon.all.loom")
relumi.all.ss3.day0 <- apply(umi.all.ss3.day0,2,umi2rel)
relumi.melt.all.ss3.day0 <- as.data.table(melt(relumi.all.ss3.day0))
relumi.melt.all.ss3.day0 %<>% add.info.ss3(coldata = meta.ss3.day0)

umi.all.ss3.cl5diff <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_diff_clone5.umicount.exon.all.loom")
relumi.all.ss3.cl5diff <- apply(umi.all.ss3.cl5diff,2,umi2rel)
relumi.melt.all.ss3.cl5diff <- as.data.table(melt(relumi.all.ss3.cl5diff))
relumi.melt.all.ss3.cl5diff %<>% add.info.ss3(coldata = meta.ss3.cl5diff)

relumi.melt.all.ss3 <- rbindlist(list(
  relumi.melt.all.ss3.day0,
  relumi.melt.all.ss3.cl5diff
))
relumi.melt.all.ss3 %<>% .[exclude == F]
relumi.melt.all.ss3[,expressed := mean(value,na.rm = T)>0, by = c("gene","day")]

# load allelic UMI
umi.melt.ss3.day0 <- load.umi.ss3("smartseq3/results/zUMIs_output/allelic/mESC_day0.BL6_direct_UMIs.txt","smartseq3/results/zUMIs_output/allelic/mESC_day0.CAST_direct_UMIs.txt", coldata=meta.ss3.day0)
umi.melt.ss3.cl5diff <- load.umi.ss3("smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.BL6_direct_UMIs.txt", "smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.CAST_direct_UMIs.txt", coldata=meta.ss3.cl5diff)

umi.melt.ss3 <- rbindlist(list(umi.melt.ss3.day0, umi.melt.ss3.cl5diff))
umi.melt.ss3 %<>% .[exclude == F]
umi.melt.ss3[,expressed := relumi.melt.all.ss3$expressed[match(paste0(gene,day),with(relumi.melt.all.ss3,paste0(gene,day)))]]

## cleanup
rm(list = c(
  "altratio.ss3.day0",
  "altratio.ss3.day0.filt",
  "altratio.ss3.cl5diff",
  "altratio.ss3.cl5diff.filt",
  "rpkm.all.ss3.day0",
  "rpkm.all.ss3.cl5diff",
  "umi.melt.ss3.day0",
  "umi.melt.ss3.cl5diff",
  "umi.all.ss3.day0",
  "umi.all.ss3.cl5diff",
  "relumi.all.ss3.day0",
  "relumi.all.ss3.cl5diff",
  "relumi.melt.all.ss3.day0",
  "relumi.melt.all.ss3.cl5diff",
  "tpm.all.ss3.day0",
  "tpm.all.ss3.cl5diff",
  "tpm.allele.ss3.day0",
  "tpm.allele.ss3.cl5diff",
  "tpm.melt.all.ss3.day0",
  "tpm.melt.all.ss3.cl5diff",
  "meta.ss3.day0",
  "meta.ss3.cl5diff",
  "counts.ss3.day0",
  "counts.ss3.cl5diff",
  "tpm.melt.ss3.day0",
  "tpm.melt.ss3.cl5diff"
))