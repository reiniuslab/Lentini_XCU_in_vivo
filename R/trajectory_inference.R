suppressMessages(source("R/load_data_cellrep.R"))

### set up functions

get.sce <- function(x, coldata, cutoff = 1){
  require(magrittr)
  require(scater)
  require(scran)
  
  # convert to SCE format
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(x)
    ),
    colData = coldata
  )
  
  # filter for expressed genes
  sce.filt <- subset(sce,calcAverage(sce) >= cutoff)
  sce.filt %<>% computeSumFactors()
  sce.filt %<>% normalize()
  
  sce.filt
}

trajectory.inference <- function(x){
  require(slingshot)
  require(mclust)
  require(magrittr)
  x %<>% runDiffusionMap()
  cl <- Mclust(reducedDim(x, "DiffusionMap"))$classification
  x %<>% slingshot(clusterLabels = cl, reducedDim = "DiffusionMap")
  x$slingPseudotime_all <- rowMeans(as.data.frame(colData(x)[,grepl("slingPseudotime",colnames(colData(x)))]), na.rm = T)
  x
}

geom_slingshot <- function(data, ... ){
  require(slingshot)
  require(ggplot2)
  require(data.table)
  scs <- slingCurves(data)
  pc <- rbindlist(lapply(slingCurves(data),function(x) with(x,data.table(s[ord,1:2])) ),idcol = T)
  return(geom_path(data = pc, aes(DC1, DC2, group=.id), ... ))
}

###

library(scater)
library(scran)
library(MAST)
# remove missing and chrY genes
idx.gn <- with(refseq.gene, !is.na(name2) & chrom != "chrY")

# convert to SCE format
sce.all <- get.sce(counts.c57 + counts.cast, meta)

sce.split.lineage <- list(
  EPI = subset(sce.all, , Lineage == "EPI"),
  ExE = subset(sce.all, , Lineage == "VE"),
  VE = subset(sce.all, , Lineage == "ExE")
)
sce.split.lineage %<>% lapply(trajectory.inference)

## write results
# fwrite(rbindlist(lapply(sce.split.lineage,function(x) with(colData(x), list(LibraryName=LibraryName,Pseudotime=slingPseudotime_all)))),"data/pseudotime.tsv",quote = F,sep = "\t")

## plot trajectories
library(ggplot2)
library(cowplot)
# Diffusionmap
ggsave2("plots/pseudotime_diffusionmap.pdf",width = 12, height = 4,
  multiplot(cols = 3,
    plotDiffusionMap(sce.split.lineage$EPI, colour_by="slingPseudotime_all", shape_by = "EmbryonicDay") + geom_slingshot(sce.split.lineage$EPI, size = 1, alpha=0.33) + labs(x="DC1",y="DC2", title="EPI") + scale_color_viridis_c(name="Pseudotime") + theme(legend.position = "top"),
    plotDiffusionMap(sce.split.lineage$ExE, colour_by="slingPseudotime_all", shape_by = "EmbryonicDay") + geom_slingshot(sce.split.lineage$ExE, size = 1, alpha=0.33) + labs(x="DC1",y="DC2", title="ExE") + scale_color_viridis_c(name="Pseudotime") + theme(legend.position = "top"),
    plotDiffusionMap(sce.split.lineage$VE, colour_by="slingPseudotime_all", shape_by = "EmbryonicDay") + geom_slingshot(sce.split.lineage$VE, size = 1, alpha=0.33) + labs(x="DC1",y="DC2", title="VE") + scale_color_viridis_c(name="Pseudotime") + theme(legend.position = "top")
  )
)

## use pseudotime for timeseries
tpm.melt[as.data.table(colData(sce.split.lineage$EPI)), pseudotime := slingPseudotime_all, on = "sample==LibraryName"]
tpm.melt[as.data.table(colData(sce.split.lineage$ExE)), pseudotime := slingPseudotime_all, on = "sample==LibraryName"]
tpm.melt[as.data.table(colData(sce.split.lineage$VE)), pseudotime := slingPseudotime_all, on = "sample==LibraryName"]

tpm.melt.avg <- tpm.melt[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","x.status.maternal","chrx","maternal","pseudotime")]

# plot
library(cowplot)
library(ggplot2)
## timeseries over pseudotime
p.xexprs.pt.epi <-
  ggplot(tpm.melt.avg[lineage == "EPI"],aes(y=V1, x=pseudotime, col=paste(chrx, l1 == tolower(maternal)), fill=paste(chrx, l1 == tolower(maternal)) )) +
    geom_point(stroke=NA) +
    coord_cartesian(ylim=c(0,100)) +
    facet_wrap(~ lineage + sex) +
    scale_color_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    scale_fill_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    labs(x = "Pseudotime", y = "Expression (TPM)") + 
    guides(fill=F) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

p.xexprs.pt.exe <-
  ggplot(tpm.melt.avg[lineage == "ExE"],aes(y=V1, x=pseudotime, col=paste(chrx, l1 == tolower(maternal)), fill=paste(chrx, l1 == tolower(maternal)) )) +
    geom_point(stroke=NA) +
    coord_cartesian(ylim=c(0,100)) +
    facet_wrap(~ lineage + sex) +
    scale_color_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    scale_fill_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    labs(x = "Pseudotime", y = "Expression (TPM)") + 
    guides(fill=F) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

p.xexprs.pt.ve <-
  ggplot(tpm.melt.avg[lineage == "VE"],aes(y=V1, x=pseudotime, col=paste(chrx, l1 == tolower(maternal)), fill=paste(chrx, l1 == tolower(maternal)) )) +
    geom_point(stroke=NA) +
    coord_cartesian(ylim=c(0,100)) +
    facet_wrap(~ lineage + sex) +
    scale_color_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    scale_fill_brewer(palette="Paired", name = NULL,labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal") ) +
    labs(x = "Pseudotime", y = "Expression (TPM)") + 
    guides(fill=F) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())
  
ggsave2("plots/pseudotime_x_expression_timeseries.pdf",width = 8, height = 4,
  plot_grid(nrow = 1,
    p.xexprs.pt.epi,
    p.xexprs.pt.exe,
    p.xexprs.pt.ve
  )
)

# linear regression on pseudotime

fit.pseudo <- 
rbindlist(fill = T, list(
  tpm.melt.avg[lineage == "EPI" & chrx == T,as.list(summary(lm(V1 ~ pseudotime + (l1 == tolower(maternal)) ))$coefficients[-1,4]), by = c("lineage","sex","x.status.maternal")],
  tpm.melt.avg[lineage != "EPI" & chrx == T,as.list(summary(lm(V1 ~ pseudotime + (l1 == tolower(maternal)) ))$coefficients[-1,4]), by = c("lineage","sex")]
  )
)
fit.pseudo.melt <- melt(fit.pseudo, measure.vars = c("pseudotime","l1 == tolower(maternal)TRUE"))

## regression coefficients
p.lm.pt <- 
  ggplot(fit.pseudo.melt,aes(y=-log10(value), x=paste(sex,x.status.maternal), fill= variable)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = -log10(1e-3), lty=3) +
  #geom_vline(xintercept = seq(1.5,5.5),lty=2,col=pal.t10[8]) +
  facet_grid(~lineage, space = "free_x", scales="free_x") +
  labs(y="-log10(P-value)") +
  scale_fill_brewer(palette="Set1", name="Allele (XmXp)", labels=c("Pseudotime","Allele")) +
  theme_cowplot() +
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/pseudotime_x_expression_lmcoef.pdf",width = 6, height = 3,
  p.lm.pt
)