suppressMessages(source("R/load_data_science.R"))

library(scater)
library(scran)

sce.blast <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts.c57.blast+counts.cast.blast),
    tpm = tpm.all.blast
  ),
  colData = meta.blast
)

sce.blast.filt <- subset(sce.blast,calcAverage(sce.blast) > 0)
sce.blast.filt %<>% computeSumFactors()
sce.blast.filt %<>% normalize()

# identify hvgs
var.fit.blast <- trendVar(sce.blast.filt, use.spikes=F)
var.decomp.blast <- decomposeVar(sce.blast.filt, var.fit.blast)
var.decomp.blast$name <- rownames(var.decomp.blast)
var.decomp.blast <- var.decomp.blast[with(var.decomp.blast,order(-bio,FDR)),]
hvgs.blast <- var.decomp.blast[var.decomp.blast$FDR < 0.05,'name']

# UMAP
set.seed(33) # reproducibility
sce.blast.filt %<>% runUMAP(feature_set = head(hvgs.blast,1000))

# graph-based clustering
sce.blast.filt$cluster <- quickCluster(sce.blast.filt,use.ranks=T,min.size=1,subset.row=head(hvgs.blast,1000),graph.fun= igraph::cluster_louvain)

## plot UMAP
library(cowplot)
p.umap.blast.lineage <- plotUMAP(sce.blast.filt, colour_by = "Lineage") + labs(x="UMAP1",y="UMAP2") + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))
p.umap.blast.cluster <- plotUMAP(sce.blast.filt, colour_by = "cluster") + stat_ellipse(show.legend = F, lty=2, aes(col=colour_by)) + labs(x="UMAP1",y="UMAP2") + scale_color_manual(values = pal.t10) + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))
p.umap.blast.cdx2 <- plotUMAP(sce.blast.filt, colour_by = "Cdx2", by_exprs_values = "tpm")+ stat_ellipse(show.legend = F, lty=2, aes(col=sce.blast.filt$cluster)) + labs(x="UMAP1",y="UMAP2") + scale_color_manual(values = pal.t10) + scale_fill_distiller(palette="Reds", name="Cdx2", direction = 1) + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))
p.umap.blast.krt18 <- plotUMAP(sce.blast.filt, colour_by = "Krt18", by_exprs_values = "tpm")+ stat_ellipse(show.legend = F, lty=2, aes(col=sce.blast.filt$cluster)) + labs(x="UMAP1",y="UMAP2") + scale_color_manual(values = pal.t10) + scale_fill_distiller(palette="Reds", name="Krt18", direction = 1) + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))
p.umap.blast.nanog <- plotUMAP(sce.blast.filt, colour_by = "Nanog", by_exprs_values = "tpm")+ stat_ellipse(show.legend = F, lty=2, aes(col=sce.blast.filt$cluster)) + labs(x="UMAP1",y="UMAP2") + scale_color_manual(values = pal.t10) + scale_fill_distiller(palette="Blues", name="Nanog", direction = 1) + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))
p.umap.blast.sox2 <- plotUMAP(sce.blast.filt, colour_by = "Sox2", by_exprs_values = "tpm")+ stat_ellipse(show.legend = F, lty=2, aes(col=sce.blast.filt$cluster)) + labs(x="UMAP1",y="UMAP2") + scale_color_manual(values = pal.t10) + scale_fill_distiller(palette="Blues", name="Sox2", direction = 1) + theme(legend.position = "top") + coord_cartesian(xlim=c(-5,5),ylim=c(-3,3))

ggsave2("plots/blast_subgrouping_unsupervised_umap.pdf", width = 4, height = 6,
multiplot(layout = matrix(1:6,ncol = 2,byrow = T),
  p.umap.blast.lineage,
  p.umap.blast.cluster,
  p.umap.blast.cdx2,
  p.umap.blast.krt18,
  p.umap.blast.nanog,
  p.umap.blast.sox2
)
)

## plot heatmap of markers
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
markers.blast <- c("Cdx2","Eomes","Gata3","Krt18","Fgfr2","Nanog","Pou5f1","Sox2","Fgf4","Klf4","Etv5","Gata4","Gata6","Sox17","Bmp4","Pdgfra","Col4a1","Sparc")
markers.z <- t(apply(tpm.all.blast[markers.blast,],1,scale))

## subcluster cluster4
sce.blast.filt$cluster2 <- NA

set.seed(42) # reproducible cluster order
sce.blast.filt$cluster2[sce.blast.filt$cluster == 4] <- 4 + kmeans(t(markers.z[,sce.blast.filt$cluster == 4]),2)$cluster

# plot
get.pal <- function(pal, f){
  f <- factor(f)
  setNames(pal[1:nlevels(f)],levels(f))
}

p.blast.heat.markers <- Heatmap(
  markers.z,
  name="Z-score",
  column_split = sce.blast.filt$cluster,
  column_title = NULL,
  cluster_column_slices = F,
  col = colorRamp2(seq(-2,2,length.out = 11),rev(brewer.pal(11,"RdBu"))),
  border = "black",
  heatmap_legend_param = list(direction = "horizontal"),
  top_annotation = HeatmapAnnotation(cluster=factor(sce.blast.filt$cluster), cluster2=factor(sce.blast.filt$cluster2), col=list(cluster=get.pal(pal.t10,sce.blast.filt$cluster),cluster2=get.pal(pal.t10,sce.blast.filt$cluster2)), annotation_legend_param = list(cluster = list(direction="horizontal", nrow=1), cluster2 = list(direction="horizontal", nrow=1)))
)

pdf("plots/blast_subgrouping_heatmap_NEW.pdf",width = 8,height = 6)
  draw(p.blast.heat.markers, merge_legend=T, heatmap_legend_side="top")
dev.off()

## write clusters
# write.table(colData(sce.blast.filt)[,c("cluster","cluster2")], "data/blast_subcluster.tsv", quote = F, sep = "\t", col.names = NA)