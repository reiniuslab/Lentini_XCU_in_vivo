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

read.list <- function(path, pattern, regex, FUN, ...){
  fls <- list.files(path, pattern, full.names = T)
  message(paste0(basename(fls), collapse = "\n"))
  ls <- lapply(fls, FUN, ... )
  names(ls) <- gsub(regex, "", fls)
  return(ls)
}

detect.outliers <- function(x, nmad=3, bound=c("lower","upper") ){
  bound <- match.arg(bound, c("lower","upper"))
  if(bound == "lower"){
    x < median(x) - mad(x)*nmad
  }else{
    x > median(x) + mad(x)*nmad
  }
}

fpkm2tpm <- function(fpkm){
  (fpkm / sum(fpkm,na.rm = T)) * 1e6
}

umi2rel <- function(umi){
  (umi / sum(umi,na.rm = T)) * 1e4
}

load.umi.ss3 <- function(x,y, ... ){
  require(data.table)
  ls <- list(
    "c57" = as.matrix(read.delim(x,row.names = 1)),
    "cast" = as.matrix(read.delim(y,row.names = 1))
  )
  dat.melt <- as.data.table(melt(ls))
  dat.melt %<>% add.info.ss3( ... )
  dat.melt[, umi_rel := umi2rel(value), by="var2"] # median(colSums(Reduce("+",ls),na.rm = T))
  return(dat.melt)
}

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

## load annotations
library(rtracklayer)
gtf <- import("smartseq3/Mus_musculus.GRCm38.97.chr.gtf.gz")
gene.anno.ss3 <- unique(with(gtf, data.table(gene_id, gene_name, chrom = as.character(seqnames), start)))
gene.anno.ss3 %<>% .[chrom %in% c(1:19,"X","Y")]

## load metadata
meta.dilutions <- fread("smartseq3/metadata_dilutions.tsv")

## load count data
counts.all.dilutions <- loom2matrix("smartseq3/results_dilutions/zUMIs_output/expression/dilutions.readcount.exon.all.loom")
rpkm.all.dilutions <- loom2matrix("smartseq3/results_dilutions/zUMIs_output/expression/dilutions.rpkm.exon.all.loom")
tpm.all.dilutions <- apply(rpkm.all.dilutions, 2, fpkm2tpm)

# detect outliers based on readcounts
meta.dilutions[, exclude := detect.outliers(log10(colSums(counts.all.dilutions)))[sample_bc]]

# load subsampled allelic ratios
altratio.dilutions <- read.list("smartseq3/results_dilutions/zUMIs_output/allelic", "dilutions.fract_CAST_reads", ".*_|.txt", read.delim, row.names=1)
names(altratio.dilutions)[1] <- "all"
altratio.dilutions.filt <- lapply(altratio.dilutions, function(x) x[which(apply(x,1,var, na.rm=T)>0),] )

# scale data by allele counts
tpm.allele.dilutions <- list(
  "c57" = as.matrix(tpm.all.dilutions[intersect(row.names(altratio.dilutions.filt$all), row.names(tpm.all.dilutions)), colnames(altratio.dilutions.filt$all)] * (1 - altratio.dilutions.filt$all)),
  "cast" = as.matrix(tpm.all.dilutions[intersect(row.names(altratio.dilutions.filt$all), row.names(tpm.all.dilutions)), colnames(altratio.dilutions.filt$all)] * altratio.dilutions.filt$all)
)

## load UMI data
umi.all.dilutions <- loom2matrix("smartseq3/results_dilutions/zUMIs_output/expression/dilutions.umicount.exon.all.loom")
relumi.all.dilutions <- apply(umi.all.dilutions, 2, umi2rel)

altratio.umi.dilutions <- read.list("smartseq3/results_dilutions/zUMIs_output/allelic", "dilutions.fract_CAST_direct_UMIs", ".*_|.txt", read.delim, row.names=1)
names(altratio.umi.dilutions)[1] <- "all"
altratio.umi.dilutions.filt <- lapply(altratio.umi.dilutions, function(x) x[which(apply(x,1,var, na.rm=T)>0),] )
  
# scale data
relumi.allele.dilutions <- list(
  "c57" = as.matrix(relumi.all.dilutions[intersect(row.names(altratio.umi.dilutions.filt$all), row.names(relumi.all.dilutions)), colnames(altratio.umi.dilutions.filt$all)] * (1 - altratio.umi.dilutions.filt$all)),
  "cast" = as.matrix(relumi.all.dilutions[intersect(row.names(altratio.umi.dilutions.filt$all), row.names(relumi.all.dilutions)), colnames(altratio.umi.dilutions.filt$all)] * altratio.umi.dilutions.filt$all)
)

# melt data
tpm.melt.dilutions <- as.data.table(melt(tpm.allele.dilutions))
tpm.melt.dilutions[meta.dilutions, c("pct.c57", "exclude") := list(pct_c57, exclude), on = "Var2 == sample_bc"]
tpm.melt.dilutions[gene.anno.ss3, chrom := chrom, on = "Var1 == gene_id"]
tpm.melt.dilutions[,expressed := mean(value,na.rm = T)>0, by = c("Var1")]

relumi.melt.dilutions <- as.data.table(melt(relumi.allele.dilutions))
relumi.melt.dilutions[meta.dilutions, c("pct.c57", "exclude") := list(pct_c57, exclude), on = "Var2 == sample_bc"]
relumi.melt.dilutions[gene.anno.ss3, chrom := chrom, on = "Var1 == gene_id"]
relumi.melt.dilutions[,expressed := mean(value,na.rm = T)>0, by = c("Var1")]

## combine data

altratio.melt.dilutions <- rbindlist(list(
  "reads" = as.data.table(melt(lapply(altratio.dilutions.filt, function(x) x[row.names(x) %in% gene.anno.ss3[chrom %in% 1:19, unique(gene_id)],] ))),
  "umi" = as.data.table(melt(lapply(altratio.umi.dilutions.filt, function(x) x[row.names(x) %in% gene.anno.ss3[chrom %in% 1:19, unique(gene_id)],] )))
), idcol = T)
altratio.melt.dilutions[meta.dilutions, c("pct.c57", "exclude") := list(pct_c57, exclude), on = "variable == sample_bc"]

refratio.dilutions.avg <- altratio.melt.dilutions[exclude == F, list(ratio.c57 = mean(1-value, na.rm=T)), by = c(".id", "variable", "L1", "pct.c57")]
refratio.dilutions.avg[, L1 := factor(L1, levels=c("10000", "20000", "50000", "100000", "200000", "500000", "1000000", "all"))]

refratio.dilutions.lm <- refratio.dilutions.avg[, with(summary(lm(ratio.c57*100 ~ pct.c57)), list("r2" = r.squared, "p" = coefficients[2,4], "sigma" = sigma)), by=c(".id", "L1")]

dat.melt.dilutions.avg <- rbindlist(list(
  "tpm" = tpm.melt.dilutions[expressed == T & chrom %in% 1:19 & exclude == F, base::mean(value,na.rm=T),by=c("Var2", "L1", "pct.c57")],
  "umi" = relumi.melt.dilutions[expressed == T & chrom %in% 1:19 & exclude == F, base::mean(value,na.rm=T),by=c("Var2", "L1", "pct.c57")]
), idcol = T)

## downsample genes

dat.melt.dilutions.all <- as.data.table(melt(as.matrix(altratio.dilutions.filt$all)))
dat.melt.dilutions.all[, chrom := gene.anno.ss3[match(dat.melt.dilutions.all$Var1, gene_id), chrom]]
dat.melt.dilutions.all[meta.dilutions, pct.c57 := pct_c57, on = "Var2 == sample_bc"]

set.seed(42)
dat.melt.dilutions.all.gnsample <- dat.melt.dilutions.all[chrom %in% 1:19, lapply(round(exp(seq(log(10),log(9e3), length.out = 50))), function(i) 1-mean(sample(value, size = i), na.rm=T)), by=c("Var2", "pct.c57")]
colnames(dat.melt.dilutions.all.gnsample)[-1:-2] <- round(exp(seq(log(10),log(9e3), length.out = 50)))

dat.melt.dilutions.all.gnsample2 <- melt(dat.melt.dilutions.all.gnsample, id.vars = c("Var2", "pct.c57"))
dat.melt.dilutions.all.gnsample2[, exclude := Var2 %in% meta.dilutions[exclude == T, sample_bc]]

## PLOT
library(ggplot2)
library(cowplot)
p.dilutions.refratio <-
ggplot(refratio.dilutions.avg, aes(x=pct.c57, y=ratio.c57*100)) +
  geom_point(stroke=NA, alpha=0.5) +
  geom_smooth(method="lm", formula = y ~ x + 0 ) +
  facet_grid(.id ~ L1) +
  geom_text(data = refratio.dilutions.lm[, paste0("R2 = ", round(r2,4), "\nP = ", signif(p, 3), "\nRSD = ", round(sigma,2) ), by=c(".id", "L1")], hjust=0, aes(x=0,y=80, label=V1)) +
  labs(x="Expected %C57", y="Observed %C57") +
  coord_cartesian(ylim=c(0,100), xlim=c(0,100)) +
  theme_cowplot() +
  theme(aspect.ratio = 1, strip.background = element_blank(), panel.border = element_rect(color="black"))

p.dilutions.lm <- 
ggplot(refratio.dilutions.lm, aes(x=L1, y=sigma^2, col=.id)) + 
  geom_line(aes(group=.id)) + 
  geom_point() +
  labs(x="Reads", y="Sigma^2") +
  scale_color_brewer(palette="Set1", name = NULL, labels=c("Reads", "UMIs")) +
  theme_cowplot() +
  theme(legend.position = "top", aspect.ratio = 1)

p.dilutions.exprs <-
ggplot(dat.melt.dilutions.avg, aes(x=pct.c57, y=V1, col=L1)) + 
  geom_point(stroke=NA) + 
  facet_grid(.id ~ ., scales = "free_y") +
  labs(x="Expected %C57", y="Expression", col=NULL) +
  scale_color_brewer(palette="Paired") +
  theme_cowplot() +
  theme(aspect.ratio = 1, strip.background = element_blank(), panel.border = element_rect(color="black"))

ggsave2("plots/dilutions_allelic.pdf", width = 20, height = 6,
plot_grid(align = "hv", axis = "trbl", nrow=1, rel_widths = c(1, 0.4, 0.3),
  p.dilutions.refratio, p.dilutions.lm, p.dilutions.exprs
))


##

p.dilutions.ngenes.frac <- 
ggplot(dat.melt.dilutions.all.gnsample2[exclude == F], aes(x=as.numeric(as.character(variable)), y=value*100, col=pct.c57)) + 
  geom_line(alpha=0.75,aes(group=Var2)) + 
  scale_color_distiller(palette="Spectral") + 
  labs(x="ngenes", y="observed %C57", col="expected %C57") + 
  scale_x_log10() + 
  theme_cowplot() + 
  theme(panel.background = element_rect(fill="black"))

p.dilutions.ngenes.rss <-
ggplot(dat.melt.dilutions.all.gnsample2[exclude == F, summary(aov(value*100 ~ pct.c57))[[1]][2,2], by="variable"], aes(x=as.numeric(as.character(variable)), y=V1)) + 
  geom_point() +
  labs(x="ngenes", y="RSS") + 
  scale_x_log10() + 
  theme_cowplot()

ggsave2("plots/dilutions_ngenes.pdf", width = 12, height = 6,
plot_grid(rel_widths = c(1,0.75),
  p.dilutions.ngenes.frac,
  p.dilutions.ngenes.rss
))
