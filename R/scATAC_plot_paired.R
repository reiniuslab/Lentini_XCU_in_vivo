suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

## 

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

# Peak to gene correlations
ap.total.filt2 %<>% addPeak2GeneLinks(reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

ls.p2g <- plotPeak2GeneHeatmap(ap.total.filt2, k = 25, groupBy = "Day", returnMatrices = T)

# Calculate average scATAC matrices
se.gene.total <- getMatrixFromProject(ap.total.filt, "GeneScoreMatrix", logFile = NULL)
se.gene.total %<>% .[,order(colnames(se.gene.total))]

se.gene.refratio <- getAllelicRatio(ap.allele.filt, "GeneCountMatrix")
se.gene.refratio$Group <- gsub("_C57|_CAST","",se.gene.refratio$Group)

se.gene.average <- SummarizedExperiment(
  assays=list(
    GeneScore = Matrix(getGroupedMatrix(assay(se.gene.total), se.gene.total$Group, rowMeans, na.rm=T), sparse = T),
    RefRatio = Matrix(getGroupedMatrix(assay(se.gene.refratio), se.gene.refratio$Group, rowMeans, na.rm=T), sparse = T),
    AllelicDetection = Matrix(getGroupedMatrix(!is.na(assay(se.gene.refratio)), se.gene.refratio$Group, rowMeans, na.rm=T), sparse = T)
  ),
  rowData = rowData(se.gene.total)
)

# melt and add metadata
se.gene.melt <- as.data.table(melt(list(
  "c57" = as.matrix(assay(se.gene.average, "GeneScore") * assay(se.gene.average, "RefRatio")),
  "cast" = as.matrix(assay(se.gene.average, "GeneScore") * (1 - assay(se.gene.average, "RefRatio")))
)))
colnames(se.gene.melt) <- c("rn", "group", "genescore", "allele")
se.gene.melt[, c("gene_name", "chrom") := list(as.character(rowData(se.gene.average)$name[rn]), as.character(rowData(se.gene.average)$seqnames)[rn])]
se.gene.melt[, chrom := factor(chrom, levels=paste0("chr", c(1:19,"X")))]
se.gene.melt[, sex := substr(group, 1, 1)]
se.gene.melt[chrom %in% paste0("chr",c(1:19,"X")), chrx := chrom == "chrX"]
se.gene.melt[, detected := melt(as.matrix(assay(se.gene.average, "AllelicDetection")))$value[rn]]
se.gene.melt[, x.status := gsub(".*_","",group)]
se.gene.melt[, x_allele := ifelse(allele == "c57",substr(x.status,1,2),substr(x.status,3,4))]

# Add scRNA average
se.gene.melt[tpm.melt.paired[atac == T & exclude == F, mean(value, na.rm=T), by=c("gene_name", "group", "allele")], tpm := V1, on = c("gene_name == gene_name", "group == group", "allele == allele")]
se.gene.melt[, rel_genescore := genescore/genescore[group == "Female_XaXa"], by = c("allele", "gene_name")]
se.gene.melt[, rel_tpm := tpm/tpm[group == "Female_XaXa"], by = c("allele", "gene_name")]

se.gene.melt[x.status != "XaXa" & is.finite(rel_tpm), bin := cut(log2(rel_tpm), quantile(log2(rel_tpm), probs = seq(0,1,0.1), na.rm=T))]

## Get allelic gene counts
se.count.allele <- getMatrixFromProject(ap.allele.filt, "GeneCountMatrix", logFile = NULL)
se.count.allele %<>% .[,order(colnames(se.count.allele))]

se.count.melt <- as.data.table(melt(as.matrix(assay(se.count.allele))))
colnames(se.count.melt) <- c("rn", "sample", "count")
se.count.melt[, c("gene_name", "chrom", "pos") := list(as.character(rowData(se.gene.average)$name[rn]), as.character(rowData(se.gene.average)$seqnames)[rn], as.character(rowData(se.gene.average)$start)[rn])]
se.count.melt[as.data.table(as.data.frame(colData(se.count.allele)), keep.rownames = T), c("allele", "names", "x.status", "sex", "group") := list(Allele, Names, x.status, Sex, Group), on = c("sample == rn")]
se.count.melt[, detected := mean(count>0, na.rm=T)>0.01, by="gene_name"]

## Plot
# ATAC-RNA heatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
cha.p2g <- HeatmapAnnotation(Day = ls.p2g$RNA$colData$groupBy, border = c(Day = T), col = list(Day = setNames(brewer.pal(4, "Set1"), c("D1", "D2", "D4", "D7"))), annotation_legend_param = list(Day = list(direction = "horizontal", nrow=1)), show_annotation_name = F )

h.p2g.ord <- order(apply(ls.p2g$ATAC$matrix, 1, function(x) which.max(tapply(x, ls.p2g$ATAC$colData$groupBy, mean)) ))

h.p2g.1 <- Heatmap(ls.p2g$ATAC$matrix, row_order = h.p2g.ord, name = "Z-score (ATAC)", border = T, top_annotation = cha.p2g, column_order = order(ls.p2g$ATAC$colData$groupBy), col = colorRamp2(seq(-4,4, length.out = 11), rev(brewer.pal(11, "PRGn")) ), use_raster = T, row_title_rot = 0, row_gap = unit(0, "mm"), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, heatmap_legend_param = list(direction = "horizontal") )
h.p2g.2 <- Heatmap(ls.p2g$RNA$matrix, row_order = h.p2g.ord, name = "Z-score (RNA)", border = T, top_annotation = cha.p2g, column_order = order(ls.p2g$RNA$colData$groupBy), col = colorRamp2(seq(-4,4, length.out = 11), rev(brewer.pal(11, "RdBu")) ), use_raster = T, row_title_rot = 0, row_gap = unit(0, "mm"), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, heatmap_legend_param = list(direction = "horizontal") )

pdf("plots/scatac_heatmap_peak2gene.pdf", width = 4, height = 6)
  draw(h.p2g.1 + h.p2g.2, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# Allelic
library(ggplot2)
library(cowplot)
library(scales)
library(ggpointdensity)

p.genescore.allelic <-
ggplot(dcast(se.gene.melt[x.status %in% c("XaXa", "XaXi", "XiXa") & chrom %in% c("chr8", "chrX")], group+gene_name+chrom~allele, value.var = "genescore"), aes(x=c57, y=cast, col=chrom)) +
  geom_point(alpha=0.5) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~group) +
  labs(x="Accessibility (C57 allele)", y="Accessibility (CAST allele)") +
  coord_cartesian(xlim = c(0,1.5), ylim = c(0,1.5) ) +
  scale_color_brewer(palette="Paired", name = NULL) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"), aspect.ratio = 1)

p.genecount.heatmap.chr8 <-
ggplot(se.count.melt[chrom == "chr8" & !is.na(group) & detected == T], aes(x=factor(pos, levels=sort(unique(pos))), y=reorder(names,count,sum), fill=log2(count+1) )) +
  geom_raster() +
  facet_grid(sex+x.status~chrom+allele, scales="free", space="free") +
  scale_fill_distiller(palette="Greys", direction = 1, limits=c(0,4), name="log2 counts + 1", oob=squish) +
  theme_void() +
  theme(legend.position = "top", panel.border = element_rect(colour = "black"))
p.genecount.heatmap.chrx <-
ggplot(se.count.melt[chrom == "chrX" & !is.na(group) & detected == T], aes(x=factor(pos, levels=sort(unique(pos))), y=reorder(names,count,sum), fill=log2(count+1) )) +
  geom_raster() +
  facet_grid(sex+x.status~chrom+allele, scales="free", space="free") +
  scale_fill_distiller(palette="Greys", direction = 1, limits=c(0,4), name="log2 counts + 1", oob=squish) +
  theme_void() +
  theme(legend.position = "top", panel.border = element_rect(colour = "black"))

p.paired.score.2d <-
ggplot(se.gene.melt[group %in% c("Female_XaXa", "Female_XaXi", "Female_XiXa", "Male_XaXi") & chrx == T], aes(log2(genescore+1), log2(tpm+1), col= allele )) +
  geom_point() +
  facet_wrap(~group) +
  labs(x="log2 accessibility + 1", y="log2 expression + 1", col = "Allele") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,10)) +
  scale_color_brewer(palette = "Paired", name = "Allele") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"), aspect.ratio = 1)

p.paired.density.2d <-
ggplot(se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi")], aes(log2(rel_genescore), log2(rel_tpm) )) +
  geom_density_2d() +
  geom_vline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi") & chrx == F, median(log2(rel_genescore), na.rm=T), by=c("allele", "x_allele", "group")],lty = 1, aes(xintercept = V1)) +
  geom_hline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi") & chrx == F, median(log2(rel_tpm), na.rm=T), by=c("allele", "x_allele", "group")],lty = 1, aes(yintercept = V1)) +
  #geom_vline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi"), median(log2(rel_genescore), na.rm=T), by=c("chrom","allele", "x_allele", "group")],lty = 2, aes(xintercept = V1)) +
  #geom_hline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi"), median(log2(rel_tpm), na.rm=T), by=c("chrom","allele", "x_allele", "group")],lty = 2, aes(yintercept = V1)) +
  facet_grid(group+x_allele~chrom) +
  labs(x="log2 relative accessibility", y="log2 relative expression") +
  coord_cartesian(xlim=c(-2,2), ylim=c(-2,2)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), aspect.ratio = 1, panel.border = element_rect(colour = "black"))

p.paired.density.2d.ybox <- 
ggplot(se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi")], aes("X", log2(rel_tpm) )) +
  geom_boxplot(outlier.stroke = NA, notch = T) +
  geom_hline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi") & chrx == F, median(log2(rel_tpm), na.rm=T), by=c("allele", "x_allele", "group")],lty = 1, aes(yintercept = V1)) +
  facet_grid(group+x_allele~chrom) +
  labs(x="log2 relative accessibility", y="log2 relative expression") +
  coord_cartesian(ylim=c(-2,2)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), aspect.ratio = 1, panel.border = element_rect(colour = "black"))

p.paired.density.2d.xbox <- 
ggplot(se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi")], aes("Y", log2(rel_genescore))) +
  geom_boxplot(outlier.stroke = NA, notch = T) +
  geom_hline(data = se.gene.melt[group %in% c("Female_XaXi", "Female_XiXa", "Male_XaXi") & chrx == F, median(log2(rel_genescore), na.rm=T), by=c("allele", "x_allele", "group")],lty = 1, aes(yintercept = V1)) +
  facet_grid(group+x_allele~chrom) +
  labs(y="log2 relative accessibility", x="log2 relative expression") +
  coord_flip(ylim=c(-2,2)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), aspect.ratio = 1, panel.border = element_rect(colour = "black"))


## 

ggsave2("plots/scatac_genescore_allelic.pdf", width = 4, height = 4, p.genescore.allelic)

ggsave2("plots/scatac_paired_rna_score_2d.pdf", width = 4, height = 4, p.paired.score.2d)

ggsave2("plots/scatac_paired_rna_density_2d.pdf", width = 32, height = 36,
  plot_grid(ncol=1, align = "hv", axis = "trbl",
    p.paired.density.2d,
    p.paired.density.2d.ybox,
    p.paired.density.2d.xbox
  )
)
  
ggsave2("plots/scatac_heatmap_allelic.pdf", width = 8, height = 8, 
  plot_grid(nrow=1,
    p.genecount.heatmap.chr8,
    p.genecount.heatmap.chrx
  )
)

# Binned
se.gene.melt[chrx == T & x.status %in% c("XaXi", "XiXa") & !is.na(bin) & x_allele == "Xa" & !is.na(log2(rel_tpm)) & is.finite(log2(rel_tpm)), summary(lm(log2(rel_tpm) ~ bin))]
se.gene.melt[chrx == T & x.status %in% c("XaXi", "XiXa") & !is.na(bin) & x_allele == "Xa" & !is.na(log2(rel_genescore)) & is.finite(log2(rel_genescore)), summary(lm(log2(rel_genescore) ~ bin))]

p.paired.bin.cnt <- ggplot(se.gene.melt[chrx == T & x.status != "XaXa" & !is.na(bin) & x_allele == "Xa"], aes(x=bin, fill=allele)) + geom_bar() + labs(x="Expression bins", y="Count") + facet_grid(~sex+x.status) + coord_cartesian() + scale_fill_brewer(palette="Paired", drop=F) + theme_cowplot() + theme(strip.background = element_blank(), axis.text.x = element_blank())
p.paired.bin.rna <- ggplot(se.gene.melt[chrx == T & x.status != "XaXa" & !is.na(bin) & x_allele == "Xa"], aes(x=bin, y=log2(rel_tpm), col=allele)) + geom_jitter() + geom_hline(yintercept = 0, lty=2, col="grey") + labs(x="Expression bins", y="log2 relative expression") + facet_grid(~sex+x.status) + coord_cartesian(ylim=c(-3,3)) + scale_color_brewer(palette="Paired", drop=F) + theme_cowplot() + theme(strip.background = element_blank(), axis.text.x = element_blank())
p.paired.bin.atac <- ggplot(se.gene.melt[chrx == T & x.status != "XaXa" & !is.na(bin) & x_allele == "Xa"], aes(x=bin, y=log2(rel_genescore), col=allele)) + geom_jitter() + geom_hline(yintercept = 0, lty=2, col="grey") + labs(x="Expression bins", y="log2 relative accessibility") + facet_grid(~sex+x.status) + coord_cartesian(ylim=c(-3,3)) + scale_color_brewer(palette="Paired", drop=F) + theme_cowplot() + theme(strip.background = element_blank(), axis.text.x = element_blank())

ggsave2("plots/scatac_paired_rna_binned.pdf", width = 6, height = 6,
  plot_grid(ncol=1, align = "v", axis = "rl",
    p.paired.bin.cnt,
    p.paired.bin.rna,
    p.paired.bin.atac
))


## Allelic tracks
# tracks
library(grid)
library(gridExtra)
pal.paired <- RColorBrewer::brewer.pal(12,"Paired")

p.tracks.chr <- plotBrowserTrack(ap.allele.filt[!is.na(ap.allele.filt$Group)], plotSummary = c("bulkTrack", "featureTrack"), groupBy = "Group", tileSize = 1e5, region = ap.allele.filt@genomeAnnotation$chromSizes[1:20], pal = pal.paired)

# grid::grid.draw(p.tracks.chr[[20]])

ggsave2("plots/scatac_tracks_chr.pdf", width = 12, height = 80, limitsize = F, 
  plot_grid(ncol = 1,plotlist = p.tracks.chr)
)

# rna
tpm.melt.paired.avg <- tpm.melt.paired[chrom %in% c(5, "X") & !is.na(group) & exclude == F & expressed == T, base::mean(value, na.rm=T, trim=0.2), by=c("sample", "allele", "group", "chrom")]

p.tracks.expr <-
ggplot(tpm.melt.paired.avg, aes(x=allele, y=V1, col=allele)) + 
  geom_boxplot() + 
  facet_grid(group~chrom) + 
  coord_cartesian(ylim=c(0,60)) +
  labs(x=NULL, y="Expression (TPM)") +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/scatac_track_chr_exprs.pdf", height = 6, width = 4, p.tracks.expr)
