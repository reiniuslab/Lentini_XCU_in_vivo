library(data.table)
library(magrittr)

### set up functions

scale.by.allele <- function(x, ref, alt, labs=c("c57","cast")){
  total <- ref+alt
  refratio <- ref/total
  x.ref <- x*refratio
  x.alt <- x*(1-refratio)
  x.allele <- list(x.ref, x.alt)
  names(x.allele) <- labs
  return(x.allele)
}

fpkm2tpm <- function(fpkm){
  (fpkm / sum(fpkm,na.rm = T)) * 1e6
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

add.info <- function(x, rowdata, coldata, rowjoin = "Var1==name2", coljoin = "Var2==LibraryName"){
  require(data.table)
  require(magrittr)
  
  x %<>% as.data.table()
  rowdata %<>% as.data.table()
  coldata %<>% as.data.table()
  
  x[rowdata, chr := chrom, on = rowjoin]
  x[chr != "chrY",chrx := chr == "chrX"]
  x %<>% .[coldata, on = coljoin]
  
  idx <- match(c(gsub("==.*","",rowjoin),gsub("==.*","",coljoin)),colnames(x))
  
  colnames(x)[idx] <- c("gene", "sample")
  colnames(x) %<>% tolower()
  return(x[!is.na(gene) | !is.na(chr)])
}

substr.sort <- function(x1,x2){
  m <- matrix(c(x1,x2),ncol = 2,byrow = F)
  apply(m,1,function(y) paste0(sort(y),collapse = "") )
}

substr.flip <- function(x1,x2,index){
  m <- matrix(c(x1,x2),ncol = 2,byrow = F)
  # switch order
  m[index,] <- t(apply(m[index,],1, rev))
  # collapse to string
  apply(m,1,paste0,collapse = "")
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

### end

# define palette
pal.t10 <- c("#1c78b6", "#f07e20", "#2ba037", "#d72728", "#8e67a9", "#8d574c", "#d379b0", "#7f7f80", "#bcbd20", "#2ab8cb")

# load data
counts.c57 <- read.delim("cell_rep/counts_c57.tsv", row.names = 1)
counts.cast <- read.delim("cell_rep/counts_cast.tsv", row.names = 1)
rpkm.all <- as.matrix(read.delim("cell_rep/rpkm_all.tsv", row.names = 1))

meta <- fread("cell_rep/metadata.tsv")
meta[,Maternal := gsub(" .*","",GeneticBackground)]

refseq.gene <- fread("cell_rep/refseq_gene_annotation.tsv")

# determine x status
c57.x.frac <- get.x.frac.ref(counts.c57, counts.cast, refseq.gene$chrom)
x.status <- get.x.status(c57.x.frac)

meta$c57.x.frac <- c57.x.frac
# Xc57Xcast
meta$x.status <- x.status
# non-allelic specific x.status (e.g. for chr:Autosome ratio)
meta[,x.status.simple := substr.sort(substr(x.status,1,2),substr(x.status,3,4))]
# XmatXpat
meta[,x.status.maternal := substr.flip(substr(x.status,1,2),substr(x.status,3,4), Maternal == "CAST")]

# get all group combinations
spl.anno <- meta[,list(sex=Sex,lineage=Lineage,xstatus=x.status)]
spl.idx <- interaction(spl.anno,drop=T)

# convert to tpm
tpm.all <- apply(rpkm.all,2,fpkm2tpm)
tpm.allele <- scale.by.allele(tpm.all, counts.c57, counts.cast)
tpm.allele %<>% lapply(as.matrix)

## melt data
# tpm total
tpm.melt.all <- melt(tpm.all)
tpm.melt.all %<>% add.info(refseq.gene, meta)
tpm.melt.all[,expressed_lineage := mean(value,na.rm = T)>0, by = c("gene","lineage")]

# tpm per allele
tpm.melt <- melt(tpm.allele)
tpm.melt %<>% add.info(refseq.gene, meta)
tpm.melt[,expressed_lineage := tpm.melt.all$expressed_lineage[match(paste0(gene,lineage),with(tpm.melt.all,paste0(gene,lineage)))]]

## cleanup
rm(list = c("c57.x.frac","x.status"))
