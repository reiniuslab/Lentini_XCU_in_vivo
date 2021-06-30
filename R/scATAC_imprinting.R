suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

## imprinted genes from doi:10.1038/ng.3678
gn.imp <- c("Sgce", "Peg3","Peg12", "Plagl1", "Zrsr1", "Peg13", "Airn", "Impact", "Nckap5", "Dlx5", "Gm5422", "Grb10", "Meg3", "Rian", "Mirg", "Igf2r", "Igf2", "H19")

## maternal ratios
se.gene.refratio.a <- getAllelicRatio(ap.allele.filt, useSeqnames = paste0("chr", 1:19), "GeneCountMatrix")
se.gene.total.a <- getMatrixFromProject(ap.total.filt, useSeqnames = paste0("chr", 1:19) )

dat.imp <- data.table(gene_name = rowData(se.gene.total.a)$name, gene_id = gene.anno.ss3[match(rowData(se.gene.total.a)$name, gene_name), gene_id], refratio_atac = rowMeans(assay(se.gene.refratio.a), na.rm=T))
dat.imp[, refratio_rna := 1-rowMeans(altratio.paired.filt[match(gene_id, row.names(altratio.paired.filt)),], na.rm=T) ]
dat.imp[, expression_rna := rowMeans(tpm.all.paired[match(gene_id, row.names(tpm.all.paired)),], na.rm=T) ]
dat.imp[, detected_rna := rowMeans(tpm.all.paired[match(gene_id, row.names(tpm.all.paired)),] > 0, na.rm=T) ]
dat.imp[, imprinted := gene_name %in% gn.imp]

## PLOT
library(ggplot2)
library(cowplot)
# tracks
p.tracks.imp <- plotBrowserTrack(ap.allele.filt, groupBy = "Allele", geneSymbol = c("Meg3", "Impact"), plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 4), scCellsMax = 2000, tileSize = 50, upstream = c(1.5e4, 8e3), downstream = c(5e3, 8e3) )

ggsave2("plots/scatac_tracks_imprinted.pdf", width = 6, height = 6, plot_grid(ncol = 1, plotlist = p.tracks.imp))

# scatter
library(ggrepel)
p.scatter.imp <-
ggplot(dat.imp, aes(x=refratio_atac, y=refratio_rna)) + 
  stat_ellipse(data = dat.imp[imprinted == F & expression_rna > 1], lty = 2) +
  geom_point(data = dat.imp[imprinted == T & expression_rna > 1], aes(size = detected_rna)) + 
  geom_text_repel(data = dat.imp[imprinted == T & expression_rna > 1], aes(label=gene_name)) +
  labs(x="Maternal ratio (ATAC)", y="Maternal ratio (RNA)", size ="Fraction detected") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/scatac_scatter_imprinted.pdf", width = 5, height = 4, p.scatter.imp)