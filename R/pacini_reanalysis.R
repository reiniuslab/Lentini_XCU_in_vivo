##

get.x.status <- function(frac){
  x <- NA
  x[frac >= 0.9] <- "XaXi"
  x[frac <= 0.1] <- "XiXa"
  x[frac < 0.9 & frac >= 0.6] <- "XaXs"
  x[frac > 0.1 & frac <= 0.4] <- "XsXa"
  x[frac > 0.4 & frac < 0.6] <- "XaXa"
  x
}

##

## gene annotations
library(data.table)
gene.anno.mm10 <- unique(fread("Pacini/mm10_gene_anno.txt.gz")[,c("name2", "chrom")])

## Total umi counts
library(scater)
library(scran)
library(magrittr)
sce.pacini <- SingleCellExperiment(assays= list(counts=as.matrix(read.delim("Pacini/GSE151009_UMICountMatrix.txt.gz")) ))

sce.pacini %<>% calculateQCMetrics()

sce.pacini.filt <- subset(sce.pacini,, !isOutlier(sce.pacini$total_counts, nmads=3, type="lower", log=T) )
sce.pacini.filt %<>% computeSumFactors()
sce.pacini.filt %<>% normalizeSCE()

## Allelic

umi.allelic.pacini <- list(
  "c57" = as.matrix(read.delim("Pacini/GSE151009_B6_UMICountMatrix.txt.gz")),
  "cast" = as.matrix(read.delim("Pacini/GSE151009_Cast_UMICountMatrix.txt.gz"))
)

umi.allelic.pacini.melt <- as.data.table(melt(umi.allelic.pacini))
umi.allelic.pacini.melt[, day := substr(Var2, 2,2)]
umi.allelic.pacini.melt[gene.anno.mm10, chrom := chrom, on=c("Var1 == name2")]
umi.allelic.pacini.melt[, reference_ratio := sum(value[L1 == "c57" & chrom == "chrX"], na.rm=T) / (sum(value[L1 == "c57" & chrom == "chrX"], na.rm=T) + sum(value[L1 == "cast" & chrom == "chrX"], na.rm=T)), by = "Var2"]

#
logumi.all.pacini <- logcounts(subset(sce.pacini.filt, calcAverage(sce.pacini.filt) > 0))

logumi.all.pacini.melt <- as.data.table(melt(logumi.all.pacini))
logumi.all.pacini.melt[, day := substr(Var2, 2,2)]
logumi.all.pacini.melt[gene.anno.mm10, chrom := chrom, on=c("Var1 == name2")]
logumi.all.pacini.melt[unique(umi.allelic.pacini.melt[,list(Var2, reference_ratio)]), reference_ratio := reference_ratio, on=c("Var2 == Var2")]
logumi.all.pacini.melt[, x.status := ifelse(reference_ratio %between% c(0.4,0.6), "XaXa", "XaXi")]

#
library(matrixStats)
refratio.pacini <- umi.allelic.pacini[["c57"]] / (umi.allelic.pacini[["c57"]] + umi.allelic.pacini[["cast"]])

idx.row <- intersect(row.names(logumi.all.pacini), row.names(refratio.pacini[rowVars(refratio.pacini, na.rm = T) > 0,]))
idx.col <- intersect(colnames(logumi.all.pacini), colnames(refratio.pacini))

logumi.pacini.melt <- na.omit(as.data.table(melt(list(
  "c57" = as.matrix(logumi.all.pacini[idx.row, idx.col]) * as.matrix(refratio.pacini[idx.row, idx.col]),
  "cast" = as.matrix(logumi.all.pacini[idx.row, idx.col]) * (1-as.matrix(refratio.pacini[idx.row, idx.col]))
))))
logumi.pacini.melt[, day := substr(Var2, 2,2)]
logumi.pacini.melt[gene.anno.mm10, chrom := chrom, on=c("Var1 == name2")]
logumi.pacini.melt[unique(umi.allelic.pacini.melt[,list(Var2, reference_ratio)]), reference_ratio := reference_ratio, on=c("Var2 == Var2")]
logumi.pacini.melt[, x.status := get.x.status(reference_ratio)]
logumi.pacini.melt[chrom != "chrY", chrx := chrom == "chrX"]

# associations
logumi.pacini.melt[chrom == "chrX" & L1 == "c57", mean(value), by=c("Var2", "reference_ratio")][, summary(lm(V1~reference_ratio))]
logumi.pacini.melt[chrom == "chrX" & L1 == "cast", mean(value), by=c("Var2", "reference_ratio")][, summary(lm(V1~reference_ratio))]

## PLOT
library(ggplot2)
library(cowplot)
p.pacini.refratio <-
ggplot(unique(umi.allelic.pacini.melt[, list(Var2, reference_ratio, day)]), aes(x=reference_ratio, y=..density..)) +
  geom_histogram() +
  labs(x="C57BL6/J ratio", y="Density") +
  facet_grid(~day) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

p.pacini.xexprs <-
ggplot(logumi.all.pacini.melt[chrom == "chrX", mean(value, na.rm=T),by = c("Var2", "day", "x.status")], aes(x=day, y=V1, col=x.status)) +
  geom_boxplot(outlier.stroke = NA) +
  geom_text(data= logumi.all.pacini.melt[chrom == "chrX", mean(value, na.rm=T),by = c("Var2", "day", "x.status")][, .N, by=c("day", "x.status")], position = position_dodge(0.75), show.legend = F, aes(label=N, y=1.5)) +
  labs(x="Day", y="Log-normalized UMIs") +
  expand_limits(y=c(0, 1.5)) +
  scale_color_brewer(palette="Set1", name="X status") +
  theme_cowplot()

p.pacini.xexprs.allelic <- 
ggplot(logumi.pacini.melt[chrom == "chrX" & x.status %in% c("XaXa", "XaXi", "XiXa"), mean(value, na.rm=T),by = c("Var2", "x.status", "L1", "x.status")], aes(x=x.status, y=V1, col=L1)) +
  geom_boxplot(outlier.stroke = NA) +
  labs(x="X status", y="Log-normalized UMIs") +
  expand_limits(y=c(0, 1.5)) +
  scale_color_brewer(palette="Paired", name="Allele") +
  theme_cowplot()

p.pacini.xexprs.allelic2 <-
ggplot(logumi.pacini.melt[!is.na(chrx), mean(value), by=c("Var2", "L1", "chrx", "reference_ratio")], aes(x=reference_ratio, y=V1, col=paste(chrx,L1=="c57"), fill=paste(chrx, L1=="c57") )) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(chrx~.) +
  labs(x="C57BL6/J ratio", y="Log-normalized UMIs") +
  scale_color_brewer(palette="Paired", name="Allele", labels = c("Autosomes:CAST", "Autosomes:C57", "ChrX:CAST", "ChrX:C57"), aesthetics = c("col","fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/pacini_reanalysis.pdf", width = 5, height = 6,
  plot_grid(ncol = 1, align = "v", axis = "trbl",
    p.pacini.refratio,
    p.pacini.xexprs,
    p.pacini.xexprs.allelic
  )
)

ggsave2("plots/pacini_allelefraction_expression.pdf", width = 4, height = 3, p.pacini.xexprs.allelic2)