library(tximport)
library(data.table)
library(magrittr)
library(biomaRt)

##

density.max <- function(x, ... ){
  # calculates the peak position in a density estimate
  if(length(na.omit(x)) < 2){
    return(-Inf)
  }else{
    dw <- density(x, ... )
    dw.max <- dw$x[which.max(dw$y)]
    return(dw.max)
  }
}

##

x.escapees <- fread("x_escapees.tsv")

## gene annotations
t2g <- fread("gencode.vM22.metadata.MGI.gz", header = F)
names(t2g) <- c("transcript", "gene","id")

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
gene.anno <- as.data.table(getBM(attributes = c("mgi_symbol", "chromosome_name","start_position"), filters = "mgi_id", values = unique(t2g$id), mart = mart))
gene.anno %<>% .[chromosome_name %in% c(1:19,"X","Y")]

## load RNA-seq data
# sample annotations
meta <- fread("public_data/public_data_meta.tsv")[order(project,sample)]
meta$genotype <- relevel(factor(meta$genotype),ref="XX")

fls <- list.files(c("public_data/PRJNA130405","public_data/PRJNA354946"),"quant.sf",recursive = T,full.names = T)
names(fls) <- sapply(strsplit(fls,"/"),"[[",5)

# rna quant
txi <- tximport(fls,type="salmon",tx2gene = t2g[,1:2,with=F])

txi.melt <- as.data.table(melt(txi$abundance))
txi.melt[gene.anno, chr := paste0("chr",chromosome_name),on="Var1 == mgi_symbol"]
txi.melt[chr != "chrY", chrx := chr == "chrX"]
txi.melt %<>% .[meta,on="Var2 == sample"]
txi.melt[,mean_exprs := mean(value), by=c("Var1","project","genotype")]

## diff exprs
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, meta, ~genotype+growth+project)
dds %<>%  DESeq(test="LRT", reduced = ~growth+project)

res <- rbindlist(idcol=T,
  list(
    "XXvXO" = as.data.table(as.data.frame(results(dds, contrast = c("genotype","XX","XO"))),keep.rownames = T),
    "XXvXY" = as.data.table(as.data.frame(results(dds, contrast = c("genotype","XX","XY"))),keep.rownames = T)
  )
)

res[gene.anno, chr := paste0("chr",chromosome_name), on = "rn == mgi_symbol"]
res[chr != "chrY",chrx := chr == "chrX"]
res[gene.anno, pos := start_position, on = "rn == mgi_symbol"]

# plot
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
p.bulk.xa <- 
  ggplot(txi.melt[mean_exprs > 1, median(value[chrx == T], na.rm=T) / median(value[chrx == F], na.rm=T), by=c("genotype", "Var2")], aes(x=genotype, y=V1)) +
    geom_boxplot() +
    geom_quasirandom(width = 0.1, stroke=NA) +  
    geom_hline(yintercept = c(1,0.5), lty=c(1,2)) +
    labs(x="Genotype", y="X:Autosomal ratio") +
    expand_limits(y=c(0,1.5)) +
    theme_cowplot()

ggsave2("plots/bulk_xaratio.pdf",p.bulk.xa,width = 2,height = 4)

p.bulk.fc <- 
  ggplot(res[!is.na(chrx) & baseMean >= 100],aes(x=log2FoldChange, y=..density.., fill=chrx, col=chrx)) +
    geom_density(alpha=0.33) +
    geom_vline(data=res[chrx==F & baseMean > 100,density.max(log2FoldChange, na.rm=T),by=.id], col="grey", lty=2, aes(xintercept=V1+1)) +
    facet_grid(~.id) +
    labs(y="Density", x="log2 fold change") +
    coord_cartesian(xlim=c(-3,3)) +
    scale_color_brewer(palette="Set1", name=NULL, labels=c("Autosomes","ChrX"), aesthetics = c("col","fill")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/bulk_foldchange.pdf",p.bulk.fc,width = 4,height = 2)