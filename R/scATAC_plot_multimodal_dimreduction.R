library(ArchR)
ap.total.filt2 <- loadArchRProject("scATAC/XupregTotalIntegrated")

## add impute weights
ap.total.filt2 %<>% addImputeWeights("LSI_Combined")

library(ggplot2)
library(cowplot)
## dimensionality reduction
# overview
p.dr1 <- ggplot(getEmbedding(ap.total.filt2, "UMAP_ATAC"), aes(x = `LSI_ATAC#UMAP_Dimension_1`, y = `LSI_ATAC#UMAP_Dimension_2`, col= paste(ap.total.filt2$Day, ap.total.filt2$Sex), shape = ap.total.filt2$Sex )) + geom_point() + labs(x="UMAP1", y="UMAP2", title="scATAC") + scale_color_brewer(palette="Paired", name="Day") + guides(shape=F) + theme_cowplot() + theme(panel.border = element_rect(color="black"), aspect.ratio = 1)
p.dr2 <- ggplot(getEmbedding(ap.total.filt2, "UMAP_RNA"), aes(x = `LSI_RNA#UMAP_Dimension_1`, y = `LSI_RNA#UMAP_Dimension_2`, col= paste(ap.total.filt2$Day, ap.total.filt2$Sex), shape = ap.total.filt2$Sex )) + geom_point() + labs(x="UMAP1", y="UMAP2", title="scRNA") + scale_color_brewer(palette="Paired", name="Day") + guides(shape=F) + theme_cowplot() + theme(panel.border = element_rect(color="black"), aspect.ratio = 1)
p.dr3 <- ggplot(getEmbedding(ap.total.filt2, "UMAP_Combined"), aes(x = `LSI_Combined#UMAP_Dimension_1`, y = `LSI_Combined#UMAP_Dimension_2`, col= paste(ap.total.filt2$Day, ap.total.filt2$Sex), shape = ap.total.filt2$Sex )) + geom_point() + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_color_brewer(palette="Paired", name="Day") + guides(shape=F) + theme_cowplot() + theme(panel.border = element_rect(color="black"), aspect.ratio = 1)
p.dr4 <- ggplot(getEmbedding(ap.total.filt2[ap.total.filt2$Sex == "Female"], "UMAP_Combined"), aes(x = `LSI_Combined#UMAP_Dimension_1`, y = `LSI_Combined#UMAP_Dimension_2`, z= abs(0.5-ap.total.filt2$c57.x.frac[ap.total.filt2$Sex == "Female"])/0.5 )) + stat_summary_hex(fun = function(x) mean(x==1) ) + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_fill_viridis_c(option="E", name="%XCI") + theme_cowplot() + theme(panel.border = element_rect(color="black"), aspect.ratio = 1)
# colour by markers
p.dr5 <- plotEmbedding(ap.total.filt2, "UMAP_Combined", colorBy = "GeneScoreMatrix", name = "Klf2", log2Norm = T, rastr = F, plotAs = "points", shape = 16, size=1.5) + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_colour_viridis_c(option = "A") + guides(colour = guide_colorbar("Klf2\n(accessibility)")) + theme_cowplot() + theme(panel.border = element_rect(colour="black"), aspect.ratio = 1)
p.dr6 <- plotEmbedding(ap.total.filt2, "UMAP_Combined", colorBy = "GeneExpressionMatrix", name = "Klf2", log2Norm = T, rastr = F, plotAs = "points", shape = 16, size=1.5) + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_colour_viridis_c() + guides(colour = guide_colorbar("Klf2\n(expression)")) + theme_cowplot() + theme(panel.border = element_rect(colour="black"), aspect.ratio = 1)
p.dr7 <- plotEmbedding(ap.total.filt2, "UMAP_Combined", colorBy = "GeneScoreMatrix", name = "Aif1l", log2Norm = T, rastr = F, plotAs = "points", shape = 16, size=1.5) + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_colour_viridis_c(option = "A") + guides(colour = guide_colorbar("Aif1l\n(accessibility)")) + theme_cowplot() + theme(panel.border = element_rect(colour="black"), aspect.ratio = 1)
p.dr8 <- plotEmbedding(ap.total.filt2, "UMAP_Combined", colorBy = "GeneExpressionMatrix", name = "Aif1l", log2Norm = T, rastr = F, plotAs = "points", shape = 16, size=1.5) + labs(x="UMAP1", y="UMAP2", title="Combined") + scale_colour_viridis_c() + guides(colour = guide_colorbar("Aif1l\n(expression)")) + theme_cowplot() + theme(panel.border = element_rect(colour="black"), aspect.ratio = 1)

ggsave2("plots/scatac_dimreduction.pdf", width = 8, height = 12,
plot_grid(ncol=2, align = "hv", axis = "trbl",
  p.dr1, p.dr2,
  p.dr3, p.dr4,
  p.dr5, p.dr6,
  p.dr7, p.dr8
))