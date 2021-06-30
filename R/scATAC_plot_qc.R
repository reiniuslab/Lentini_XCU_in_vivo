suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

## 

getFilteredTSS <- function(ArchRProj, flank = 1000, probs = c(0.1, 0.9)){
  require(ArchR)
  require(data.table)
  tss <- getTSS(ArchRProj) + flank
  frags <- unlist(getFragmentsFromProject(ArchRProj, logFile = NULL))
  insertions <- resize(frags, 1, fix = "start")
  
  cnt <- countOverlaps(tss, frags)
  tss[cnt %between% quantile(cnt, probs = probs)] - flank
}

##

dat.tssfrag <- as.data.table(readRDS("scATAC/QualityControl/scATAC_Total/scATAC_Total-Pre-Filter-Metadata.rds"))
dat.fraglen <- as.data.table(plotFragmentSizes(ap.total.filt2, returnDF = T))
dat.tss <- as.data.table(plotTSSEnrichment(ap.total.filt2, TSS = getFilteredTSS(ap.total.filt2), flank = 1e3, returnDF = T))
dat.ccd <- rbind(as.data.table(getCellColData(ap.total.filt2)), as.data.table(getCellColData(ap.allele.filt)), fill=T)
dat.refratio <- data.table(refratio = colMeans(assay(getAllelicRatio(ap.allele.filt, "GeneScoreMatrix", useSeqnames="chrX")), na.rm=T)[gsub(".*#", "", ap.total.filt2$cellNames)], group = ap.total.filt2$Group, refratio_rna = ap.total.filt2$c57.x.frac )

se.gene.allele <- getMatrixFromProject(ap.allele.filt, useMatrix = "GeneScoreMatrix", useSeqnames = paste0("chr", 1:19), logFile = NULL)
se.gene.allele %<>% .[,order(colnames(se.gene.allele))]

dat.geneavg <- as.data.table(getGroupedMatrix(assay(se.gene.allele), se.gene.allele$Allele, rowMeans, na.rm=T))

## PLOT
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
p.scatac.qc.tssfrag <-
ggplot(dat.tssfrag, aes(x=log10(nFrags), y=TSSEnrichment)) +
  geom_point() +
  geom_hline(yintercept = 4, lty=2) +
  geom_vline(xintercept = log10(1e3), lty=2) +
  labs(x = "log10 unique fragments", y= "TSS enrichment", col= "Neighbors") +
  scale_color_viridis_c() +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.scatac.qc.fraglen <-
ggplot(dat.fraglen, aes(x=fragmentSize, y=fragmentPercent)) +
  geom_line(lwd=1) +
  labs(x="Fragment size (bp)", y="Percentage (%)") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.scatac.qc.tss <-
ggplot(dat.tss, aes(x=x, y=normValue)) +
  geom_line() +
  labs(x="Distance to TSS (bp)", y="Normalized insertions") +
  coord_cartesian(ylim = c(0,10)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.scatac.qc.nfrags <-
ggplot(dat.ccd[, sum(nFrags[!is.na(Allele)], na.rm=T) / nFrags[is.na(Allele)], by=Names], aes(x=V1, y = ..density..)) +
  geom_histogram(fill="black") +
  labs(x="Fraction allelic fragments", y="Density") +
  expand_limits(x=0) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.scatac.qc.autoratio <-
ggplot(dat.geneavg, aes(x=log2(C57), y=log2(CAST) )) +
  geom_hex(bins=75, aes(fill=stat(density))) +
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=-8, y=3, hjust=0, label = dat.geneavg[C57 > 0 & CAST > 0, paste("rho =",round(cor(log2(C57), log2(CAST), use = "na.or.complete", method="spearman"),2))]) +
  labs(x = "log2 accessibility score (C57)", y = "log2 accessibility score (CAST)", col="Neighbors") +
  coord_cartesian(xlim=c(-10,5), ylim=c(-10,5)) +
  scale_fill_viridis_c() +
  theme_cowplot() +
  theme(aspect.ratio = 1)

# p.scatac.qc.refratio <- ggplot(dat.refratio[!is.na(group)], aes(x=gsub("_","\n",group), y=refratio)) +  geom_violin(scale="width", lty=2) +  geom_quasirandom() +  geom_text(data=dat.refratio[!is.na(group), .N, by=group], aes(label=N, y=1)) +   labs(y="scATAC maternal ratio", x=NULL) +  scale_x_discrete(limits = c("Female\nXaXa", "Female\nXaXs", "Female\nXaXi", "Female\nXsXa", "Female\nXiXa", "Male\nXaXi")) +  theme_cowplot() +  theme(aspect.ratio = 1)
p.scatac.qc.refratio <- 
ggplot(dat.refratio[!is.na(group)], aes(x=refratio_rna, y=refratio)) +
  geom_point() +
  annotate("text", x=0.2, y=0.8, label=dat.refratio[!is.na(group), with(cor.test(refratio_rna, refratio, method = "spearman"), paste0("rho = ",round(estimate,2)))]) +
  labs(y="Maternal X-fraction (ATAC)", x="Maternal X-fraction (RNA)") + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)

ggsave2("plots/scatac_qc_plots.pdf",width = 12, height = 8,
plot_grid(align = "hv", axis = "lrtb",
  p.scatac.qc.tssfrag,
  p.scatac.qc.fraglen,
  p.scatac.qc.tss,
  p.scatac.qc.nfrags,
  p.scatac.qc.autoratio,
  p.scatac.qc.refratio
))
