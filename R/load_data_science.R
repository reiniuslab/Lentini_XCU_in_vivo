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
counts.c57.blast <- read.delim("science/ss2_counts_c57.tsv", row.names = 1)
counts.cast.blast <- read.delim("science/ss2_counts_cast.tsv", row.names = 1)
rpkm.all.blast <- as.matrix(read.delim("science/ss2_rpkm_all.tsv", row.names = 1))

meta.blast <- fread("science/ss2_metadata.tsv")
meta.blast$Lineage <- factor(meta.blast$Lineage, levels=c("earlyblast","midblast","lateblast"))
meta.blast[,Maternal := gsub(" .*","",GeneticBackground)]

refseq.gene.blast <- fread("science/ss2_refseq_gene_annotation.tsv")

## scale rpkm matrix by allelic counts
#rpkm.allele.blast <- scale.by.allele(rpkm.all.blast, counts.c57.blast, counts.cast.blast)
#rpkm.allele.blast %<>% lapply(as.matrix)

# convert to tpm
tpm.all.blast <- apply(rpkm.all.blast,2,fpkm2tpm)
tpm.allele.blast <- scale.by.allele(tpm.all.blast, counts.c57.blast, counts.cast.blast)
tpm.allele.blast %<>% lapply(as.matrix)

# melt data
# tpm total
tpm.melt.all.blast <- melt(tpm.all.blast)
tpm.melt.all.blast %<>% add.info(refseq.gene.blast, meta.blast)
tpm.melt.all.blast[,expressed_lineage := mean(value,na.rm = T)>0, by = c("gene","lineage")]

# tpm per allele
tpm.melt.blast <- melt(tpm.allele.blast)
tpm.melt.blast %<>% add.info(refseq.gene.blast, meta.blast)
tpm.melt.blast[,expressed_lineage := tpm.melt.all.blast$expressed_lineage[match(paste0(gene,lineage),with(tpm.melt.all.blast,paste0(gene,lineage)))]]
