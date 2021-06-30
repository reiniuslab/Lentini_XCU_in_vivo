suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

## functions

formatPeakAnnoEnrichment <- function(x){
  x <- data.frame(TF = rownames(x), mlog10Padj = assay(x)[,1])
  x <- x[order(x$mlog10Padj, decreasing = TRUE),]
  x$rank <- seq_len(nrow(x))
  return(x)
}

## Motif enrichment 
# add motifs
ap.total.filt2 %<>% addMotifAnnotations(motifSet = "cisbp", name = "Motif")
ap.total.filt2 %<>% addBgdPeaks()
ap.total.filt2 %<>% addDeviationsMatrix(peakAnnotation = "Motif", force = TRUE)

# motif enrichment
#se.marker.peaks <- getMarkerFeatures(ap.total.filt, useMatrix = "TileMatrix", groupBy = "Day")
#dat.motifs.up <- formatPeakAnnoEnrichment(peakAnnoEnrichment(se.marker.peaks, ap.total.filt, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5"))
#dat.motifs.dn <- formatPeakAnnoEnrichment(peakAnnoEnrichment(se.marker.peaks, ap.total.filt, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5"))

# correlation with gene expression
dat.vardev <- as.data.table(getVarDeviations(ap.total.filt2, name = "MotifMatrix", plot = F))
dat.vardev[, sd2 := combinedVars > mean_sdl(combinedVars)$ymax]

se.motifs <- getGroupSE(ap.total.filt2, useMatrix = "MotifMatrix", groupBy = "Day")
se.motifs.z <- se.motifs[rowData(se.motifs)$seqnames=="z",]
rowData(se.motifs.z)$maxDelta <- lapply(seq_len(ncol(se.motifs.z)), function(x){rowMaxs(assay(se.motifs.z) - assay(se.motifs.z)[,x])}) %>% Reduce("cbind", .) %>% rowMaxs

cor.gim.mm <- as.data.table(correlateMatrices(ArchRProj = ap.total.filt2, useMatrix1 = "GeneExpressionMatrix", useMatrix2 = "MotifMatrix", reducedDims = "LSI_Combined"))
cor.gim.mm[, maxDelta := rowData(se.motifs.z)[match(MotifMatrix_name, rowData(se.motifs.z)$name), "maxDelta"] ]
cor.gim.mm %<>% .[order(abs(cor.gim.mm$cor), decreasing = TRUE)] %>% .[which(!duplicated(gsub("\\-.*","",MotifMatrix_name)))]
cor.gim.mm[, TFRegulator := cor > 0.5 & padj < 0.01 & maxDelta > quantile(maxDelta, 0.75)]

## PLOT
library(ggplot2)
library(cowplot)
# marker features
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
se.markers.atac <- getMarkerFeatures(
  ArchRProj = ap.total.filt2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Day"
)

dat.markers.atac <- plotMarkerHeatmap(se.markers.atac, plotLog2FC = T, returnMatrix = T)

idx.gns.label.atac <- c("Dppa5a", "Klf2", "Fgf4", "Cdx4", "Osm", "Cdh2", "Esrrb", "Nanog", "Lefty1", "Neil2", "Krt18", "Krt7", "Foxa2", "Rorb", "Rorc", "Rfx5", "Pou3f1", "Sox14", "Ets1", "Myof")

ha <- rowAnnotation(foo = anno_mark(at = match(idx.gns.label.atac, row.names(dat.markers.atac)), labels = idx.gns.label.atac ))

pdf("plots/scatac_heatmap_markers.pdf", width = 3, height = 4)
  Heatmap(as.matrix(dat.markers.atac), cluster_columns = F, show_row_names = F, name = "log2FC", row_order = order(apply(dat.markers.atac, 1, which.max)), right_annotation = ha, col = colorRamp2(seq(-2,2, length.out = 11), rev(brewer.pal(11, "PRGn"))) )
dev.off()
  
# compare with SS3 data
ss3.diffexpr <- fread("data/ss3_diffexpr.tsv")
ss3.diffexpr.sig <- ss3.diffexpr[fdr < 0.01 & abs(logFC) >= 0.5, name]

mat.cont <- matrix(
  c(length(union(rowData(se.markers.atac)$name, ss3.diffexpr$name)) - length(union(row.names(dat.markers.atac), ss3.diffexpr.sig)),
    length(setdiff(row.names(dat.markers.atac), ss3.diffexpr.sig)),
    length(setdiff(ss3.diffexpr.sig, row.names(dat.markers.atac))),
    length(intersect(row.names(dat.markers.atac), ss3.diffexpr.sig))
    ), nrow=2)
fisher.test(mat.cont)

dat.markers.comb <- ss3.diffexpr[fdr < 0.01 & abs(logFC) >= 0.5 & name %in% intersect(name, row.names(dat.markers.atac)), list(gene_name = name, logFC_rna = logFC, logFC_atac = dat.markers.atac[name,"D7"])]

p.markers.paired <-
ggplot(dat.markers.comb, aes(y = logFC_rna, x = logFC_atac)) +
  geom_point(stroke=NA) +
  geom_smooth(method="lm") +
  annotate("text", x = -1.5, y = 1, hjust=0, vjust=1, label = dat.markers.comb[, with(cor.test(logFC_rna, logFC_atac, method = "spearman"), paste("Rho =", round(estimate,2), "\nP =", signif(p.value,2)))]) +
  labs(x = "logFC (scATAC markers)", y = "logFC (scRNA markers)") +
  coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
  theme_cowplot() +
  theme(panel.border = element_rect(colour="black"), aspect.ratio = 1)

ggsave2("plots/scatac_paired_markers.pdf", p.markers.paired, width = 3, height = 3)

## tracks
p.tracks.markers <- plotBrowserTrack(ap.total.filt2, groupBy = "Day", geneSymbol = c("Dppa5a", "Cdx4", "Pou3f1", "Krt18"), plotSummary = c("bulkTrack", "scTrack", "featureTrack", "geneTrack"), sizes = c(10, 10, 1.5, 4), scCellsMax = 350, tileSize = 50, upstream = 1.25e4, downstream = 8.5e3)

ggsave2("plots/scatac_tracks_markers.pdf", width = 8, height = 8, limitsize = F, 
  plot_grid(plotlist = p.tracks.markers)
)

## motif enrichment
library(ggrepel)
p.motif.vardev <-
ggplot(dat.vardev[sd2 == T], aes(x=rank, y=combinedVars, col=log(rank) )) +
  geom_point(show.legend = F) +
  geom_text_repel(data=head(dat.vardev, 20), box.padding = 0.5, nudge_x = 100, aes(label=gsub("_.*", "", name), col=NULL)) +
  labs(y="Motif variability", x="Rank") +
  scale_color_viridis_c(direction = -1) +
  theme_cowplot() +
  theme(panel.border = element_rect(colour = "black"), aspect.ratio = 1)

p.motif.cor <-
ggplot(cor.gim.mm, aes(cor, maxDelta, col=cor )) +
  geom_point(show.legend = F) + 
  geom_label_repel(data = cor.gim.mm[order(maxDelta, cor, decreasing = T)][TFRegulator == T][1:10], show.legend = F, box.padding = 0.5, nudge_x = -0.25, aes(label=GeneExpressionMatrix_name, col=NULL)) +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_distiller(palette="RdBu") +
  xlab("Correlation to gene expression") +
  ylab("Max TF motif delta") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(cor.gim.mm$maxDelta)*1.05)) +
  theme_cowplot() +
  theme(panel.border = element_rect(colour = "black"), aspect.ratio = 1)

library(ggseqlogo)
ls.motifs <- lapply(mouse_pwms_v2[cor.gim.mm[TFRegulator == T][order(maxDelta, decreasing = T)]$MotifMatrix_idx], function(x) exp(as.matrix(x)) )
names(ls.motifs) <- cor.gim.mm[TFRegulator == T][order(maxDelta, decreasing = T), MotifMatrix_matchName]

p.motif.cor.seqlogo <- ggseqlogo(ls.motifs, ncol = 2) + theme_void()

#

ggsave2("plots/scatac_motif_variability.pdf", width= 6, height = 4, p.motif.vardev)

ggsave2("plots/scatac_motif_correlation.pdf", width = 7, height = 4,
  plot_grid(rel_widths = c(1.5,1),
    p.motif.cor,
    p.motif.cor.seqlogo
  )
)
