library(data.table)
library(magrittr)

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

reformat.zlm <- function(fit, hypothesis){
  require(MAST)
  res <- summary(fit, doLRT=hypothesis)
  res.dt <- res$datatable
  hurdle <- merge(res.dt[contrast==hypothesis & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  res.dt[contrast==hypothesis & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  hurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  setnames(hurdle, c("primerid","coef","Pr(>Chisq)"), c("gene","logFC","p-value"))
  setorder(hurdle,fdr)
  return(hurdle)
}

##

# load annotations
library(rtracklayer)
gtf <- import("smartseq3/Mus_musculus.GRCm38.97.chr.gtf.gz")
gene.anno.ss3 <- unique(with(gtf, data.table(gene_id, gene_name, chrom = as.character(seqnames), start)))
gene.anno.ss3 %<>% .[chrom %in% c(1:19,"X")]

meta.ss3 <- rbindlist(list(
  fread("smartseq3/metadata_ss3_day0.tsv"),
  fread("smartseq3/metadata_ss3_clone5diff.tsv")
))

counts.all.day0.ss3 <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_day0.readcount.exon.all.loom")
counts.all.cl5diff.ss3 <- loom2matrix("smartseq3/results/zUMIs_output/expression/mESC_diff_clone5.readcount.exon.all.loom")

row.idx <- intersect(row.names(counts.all.day0.ss3),row.names(counts.all.cl5diff.ss3))

library(scater)
library(scran)
library(cowplot)
sce.ss3 <- SingleCellExperiment(
  assays=list(counts=cbind(counts.all.day0.ss3[row.idx,],counts.all.cl5diff.ss3[row.idx,])),
  colData=meta.ss3[match(c(colnames(counts.all.day0.ss3),colnames(counts.all.cl5diff.ss3)),sample_bc)],
  rowData=gene.anno.ss3[match(row.idx,gene_id),1:3,with=F]
)

sce.ss3 %<>% calculateQCMetrics()
# exclude based on 3MAD read depth
sce.ss3$exclude <- isOutlier(sce.ss3$total_counts, nmads = 3, type = "lower", log = T)
## write exclusion list
# write.table(colData(sce.ss3)[,c("sample_id","exclude")],"smartseq3/samples_exclude.tsv",quote = F,sep = "\t",row.names = F)

sce.filt.ss3 <- subset(sce.ss3,,!exclude)
sce.filt.ss3 %<>% subset(calcAverage(sce.filt.ss3) > 0)
sce.filt.ss3 %<>% computeSumFactors()
sce.filt.ss3 %<>% normalize()

assay(sce.filt.ss3,"scaled_logcounts") <- t(apply(assay(sce.filt.ss3,"logcounts"),1,scale))

# identify hvgs
var.fit.ss3 <- trendVar(sce.filt.ss3, use.spikes=F)
var.decomp.ss3 <- decomposeVar(sce.filt.ss3, var.fit.ss3)
var.decomp.ss3$gene_id <- rownames(sce.ss3)
var.decomp.ss3$gene_name <- rowData(sce.ss3)$gene_name
var.decomp.ss3 <- var.decomp.ss3[with(var.decomp.ss3,order(-bio,FDR)),]
hvgs.ss3 <- var.decomp.ss3[var.decomp.ss3$FDR < 0.05,'gene_id']
names(hvgs.ss3) <- var.decomp.ss3[var.decomp.ss3$FDR < 0.05,'gene_name']

# diffusionmap of top 1,000 HVGs
p.diffmap.ss3 <- plotDiffusionMap(sce.filt.ss3, colour_by="day", shape_by="sex",rerun=T, run_args = list(feature_set = head(hvgs.ss3,1e3))) + labs(x="DC1",y="DC2")

ggsave2("plots/ss3_diffusionmap.pdf",width = 4, height = 3, p.diffmap.ss3)

# UMAP day 1
p.umap.d0.ss3 <- plotUMAP(subset(sce.filt.ss3,,day=="0"), colour_by="day", shape_by="sex",rerun=T, run_args = list(feature_set = head(hvgs.ss3,1e3))) + labs(x="UMAP1",y="UMAP2")

ggsave2("plots/ss3_umap_d0.pdf", width = 4, height = 3, p.umap.d0.ss3)

## differential expression
library(MAST)
sca.ss3 <- SceToSingleCellAssay(sce.filt.ss3)
sca.filt.ss3 <- sca.ss3[freq(sca.ss3)>0.1]

zlm.fit.ss3 <- zlm(~day+sex, sca.filt.ss3)
zlm.hurdle.ss3 <- reformat.zlm(zlm.fit.ss3, "day")
zlm.hurdle.ss3[gene.anno.ss3, name := gene_name, on = "gene==gene_id"]

# fwrite(zlm.hurdle.ss3, "data/ss3_diffexpr.tsv",quote = F,sep = "\t")

#
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
gnls <- c("Dppa3","Dppa5a","Dusp6","Eomes","Etv4","Fgf5","Fgfr1","Foxn2","Foxp1","Kif3a","Klf4","Krt18","Lefty2","Lmo4","Ly6e","Nanog","Nr0b1","Prdm14","Rex2","Sox2","Sox11","Tcl1","Usp9x","Zfp42")

mat.ss3 <- assay(sca.filt.ss3,"scaled_logcounts")[zlm.hurdle.ss3[fdr < 1e-3,head(gene,500)],]

p.heat.diffexprs.ss3 <- 
  Heatmap(
    matrix = mat.ss3,
    name = "Z-score",
    column_split = with(colData(sca.filt.ss3),factor(paste(sex,day), levels=c("M 0",paste("F",c(0,1,2,4,7)))) ),
    border = "black",
    col = colorRamp2(seq(-2,2,length.out=11), rev(brewer.pal(11,"RdBu"))),
    right_annotation = rowAnnotation(foo = anno_mark(at = match(gnls, zlm.hurdle.ss3[fdr < 1e-3,head(name,500)]), labels = gnls)),
    use_raster = T,
    cluster_columns = F,
    show_column_names = F,
    show_row_names = F,
    show_row_dend = F
  )

pdf("plots/ss3_heatmap_diffexpr.pdf",width = 5, height = 5)
  p.heat.diffexprs.ss3
dev.off()

## GSEA
# GO mapping obtained from http://www.informatics.jax.org/downloads/reports/index.html#go
go.anno <- fread("GO/go_terms.mgi",header = F)
gsc <- fread("GO/gene_association.mgi.gz", skip=24)
gsc.bp <- gsc[V9 == "P"]

sets <- with(unique(gsc.bp[,list(V3,V5)]),split(V3,V5))
sets.idx <- limma::ids2indices(sets, mcols(sca.filt.ss3)$gene_name)
sets.idx.filt <- sets.idx[lengths(sets.idx) %between% c(15,500)] # keep sets between 15 & 500 genes

#boots <- bootVcov1(zlm.fit.ss3, 100)
#saveRDS(boots,"GO/ss3_bootstraps_100.rds")
boots <- readRDS("GO/ss3_bootstraps_100.rds")

gsea <- gseaAfterBoot(zlm.fit.ss3, boots, sets.idx.filt, CoefficientHypothesis("day")) 
z.comb <- summary(gsea)
z.comb[go.anno, term := V3, on="set==V2"]
sigModules <- z.comb[combined_adj<1e-3]

# fwrite(sigModules[,list(GO=set,Term=term,Zscore=combined_Z,FDR=combined_adj)], "GO/GO_enrichment_ss3.tsv",sep = "\t")

library(ggplot2)
library(cowplot)
p.go.ss3 <-
  ggplot(sigModules[order(combined_Z)][c(1:7,(length(combined_Z)-6):length(combined_Z))], aes(x=reorder(term,combined_Z), y=combined_Z)) +
    geom_col() +
    labs(y="mESC <-- Z-score --> EPI") +
    coord_flip(ylim = c(-25,25)) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())

ggsave("plots/ss3_go_top.pdf",width = 8,height = 3,p.go.ss3)