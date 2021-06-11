library(data.table)
library(magrittr)

# load annotations
gene.anno <- fread("cell_rep/refseq_gene_annotation.tsv")

meta.invivo <- rbindlist(fill = T,idcol = T,
  list(
    "early" = fread("science/ss1_metadata.tsv"),
    "blast" = fread("science/ss2_metadata.tsv"),
    "post" = fread("cell_rep/metadata.tsv")
  )
)

# load data
counts.all.early <- read.delim("science/ss1_counts_c57.tsv",row.names = 1)+read.delim("science/ss1_counts_cast.tsv",row.names = 1)
counts.all.blast <- read.delim("science/ss2_counts_c57.tsv",row.names = 1)+read.delim("science/ss2_counts_cast.tsv",row.names = 1)
counts.all <- read.delim("cell_rep/counts_c57.tsv",row.names = 1)+read.delim("cell_rep/counts_cast.tsv",row.names = 1)

row.idx <- Reduce(intersect, list(row.names(counts.all), row.names(counts.all.blast), row.names(counts.all.early)))

library(scater)
library(scran)
# full dataset
sce.invivo <- SingleCellExperiment(
  assays=list(counts=as.matrix(cbind(counts.all.early[row.idx,], counts.all.blast[row.idx,], counts.all[row.idx,]))),
  colData=meta.invivo
)

sce.filt.invivo <- subset(sce.invivo, calcAverage(sce.invivo) > 0)
sce.filt.invivo %<>% computeSumFactors()
sce.filt.invivo %<>% normalizeSCE()
sce.filt.invivo$Lineage2 <- factor(with(colData(sce.filt.invivo), ifelse(is.na(EmbryonicDay), Lineage, EmbryonicDay)), levels=c("MIIoocyte", "zy", "early2cell", "mid2cell", "late2cell", "4cell", "8cell", "16cell", "earlyblast", "midblast", "lateblast", "5.5", "6", "6.5"))
sce.filt.invivo$Shape2 <- factor(with(colData(sce.filt.invivo), ifelse(is.na(EmbryonicDay), "EPI", Lineage)))

# only post-implantation data
sce.post <- SingleCellExperiment(assays=list(counts=as.matrix(counts.all)), colData=meta.invivo[.id == "post"])

sce.filt.post <- subset(sce.post, calcAverage(sce.post) > 0)
sce.filt.post %<>% computeSumFactors()
sce.filt.post %<>% normalizeSCE()

## identify hvgs
# full
var.fit.invivo <- trendVar(sce.filt.invivo, use.spikes=F)
var.decomp.invivo <- decomposeVar(sce.filt.invivo, var.fit.invivo)
var.decomp.invivo$gene_name <- rownames(var.decomp.invivo)
var.decomp.invivo <- var.decomp.invivo[with(var.decomp.invivo,order(-bio,FDR)),]
hvgs.invivo <- var.decomp.invivo[var.decomp.invivo$FDR < 0.05,'gene_name']

# post-implantation
var.fit.post <- trendVar(sce.filt.post, use.spikes=F)
var.decomp.post <- decomposeVar(sce.filt.post, var.fit.post)
var.decomp.post$gene_name <- rownames(var.decomp.post)
var.decomp.post <- var.decomp.post[with(var.decomp.post,order(-bio,FDR)),]
hvgs.post <- var.decomp.post[var.decomp.post$FDR < 0.05,'gene_name']

## dim reduction
# full data: PCA of top 2,000 HVGs
sce.filt.invivo %<>% runPCA(feature_set = head(hvgs.invivo,2000))

library(cowplot)
p.pca.invivo <- plotPCA(sce.filt.invivo, colour_by="Lineage2", shape_by="Shape2") + scale_color_manual(values=scales::colour_ramp(RColorBrewer::brewer.pal(11,"Spectral"))(seq(0,1,length.out=14)))

ggsave2("plots/invivo_pca.pdf",width = 4, height = 3, p.pca.invivo)

# post-implantation only:
p.pca.post <- plotPCA(sce.filt.post, colour_by = "Lineage", shape_by = "EmbryonicDay", rerun = T, run_args = list(feature_set = head(hvgs.post, 2e3))) + stat_ellipse(aes(group=colour_by))

ggsave2("plots/invivo_pca_postimplantation.pdf", width = 4, height = 3, p.pca.post)