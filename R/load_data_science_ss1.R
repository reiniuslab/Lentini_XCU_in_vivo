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
  x.avg <- rowMeans(total[chrx,],na.rm = T)
  idx <- x.avg < quantile(x.avg[x.avg > cutoff],probs = prob)
  frac <- colSums(ref[chrx,][idx,],na.rm = T)/colSums(total[chrx,][idx,],na.rm = T)
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
  return(x[!is.na(gene)])
}

### end

# define palette
pal.t10 <- c("#1c78b6", "#f07e20", "#2ba037", "#d72728", "#8e67a9", "#8d574c", "#d379b0", "#7f7f80", "#bcbd20", "#2ab8cb")

# load data
counts.c57.early <- read.delim("science/ss1_counts_c57.tsv", row.names = 1,check.names = F) # check.names=F fixes names starting with number, e.g. 16cell
counts.cast.early <- read.delim("science/ss1_counts_cast.tsv", row.names = 1, check.names = F)
rpkm.all.early <- as.matrix(read.delim("science/ss1_rpkm_all.tsv", row.names = 1, check.names = F))

meta.early <- fread("science/ss1_metadata.tsv")
meta.early$Lineage <- factor(meta.early$Lineage, levels=unique(meta.early$Lineage))
meta.early[,Maternal := meta.early[, gsub(" .*","",GeneticBackground)]]

refseq.gene.early <- fread("science/ss1_refseq_gene_annotation.tsv")

## scale rpkm matrix by allelic counts
#rpkm.allele.early <- scale.by.allele(rpkm.all.early, counts.c57.early, counts.cast.early)
#rpkm.allele.early %<>% lapply(as.matrix)

# convert to tpm and scale by allelic counts
tpm.all.early <- apply(rpkm.all.early,2,fpkm2tpm)
tpm.allele.early <- scale.by.allele(tpm.all.early, counts.c57.early, counts.cast.early)
tpm.allele.early %<>% lapply(as.matrix)

# melt data
# tpm total
tpm.melt.all.early <- melt(tpm.all.early)
tpm.melt.all.early %<>% add.info(refseq.gene.early, meta.early)
tpm.melt.all.early[,expressed_lineage := mean(value,na.rm = T)>0, by = c("gene","lineage")]

# tpm per allele
tpm.melt.early <- melt(tpm.allele.early)
tpm.melt.early %<>% add.info(refseq.gene.early, meta.early)
tpm.melt.early[,expressed_lineage := tpm.melt.all.early$expressed_lineage[match(paste0(gene,lineage),with(tpm.melt.all.early,paste0(gene,lineage)))]]