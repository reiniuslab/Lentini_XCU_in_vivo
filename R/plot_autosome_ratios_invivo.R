suppressMessages(source("R/load_data_cellrep.R"))

### functions

ratio.boot <- function(x, chr, cutoff = 1, n = 1e3){
  # column-wise
  apply(x,2,function(y){
    # filter non-expressed genes
    chr.filt <- chr[y > cutoff]
    y.filt <- y[y > cutoff]
    idx.chrx <- chr.filt == "chrX"
    idx.chry <- chr.filt == "chrY"
    idx.auto <- !chr.filt %in% c("chrX","chrY") & !is.na(chr.filt)
    median.x <- median(y.filt[idx.chrx],na.rm = T)
    median.y <- median(y.filt[idx.chry],na.rm = T)
    median.auto <- median(y.filt[idx.auto],na.rm = T)
    ratio.x <- replicate(n, median.x / median(sample(y.filt[idx.auto],size = length(median.x),replace = T)) )
    ratio.y <- replicate(n, median.y / median(sample(y.filt[idx.auto],size = length(median.y),replace = T)) )
    ratio.auto <- replicate(n, median.auto / median(sample(y.filt[idx.auto],size = length(median.auto),replace = T)) )
    return(c("Autosome"=median(ratio.auto),"chrX"=median(ratio.x), "chrY"=median(ratio.y)))
  })
}

###

## check chrX expression distribution
library(ggplot2)
library(cowplot)
p.x.ecfd <-
  ggplot(tpm.melt.all[x.status == "XaXa",mean(value), by=c("gene","chr")][!chr %in% c("chrY",NA)], aes(x=log(V1+1), col=chr=="chrX", group=chr)) +
    stat_ecdf() +
    scale_color_manual(name="EPI XaXa cells: chrX", values = c("black","red")) +
    labs(y="Ecdf", x="Expression (log(TPM+1))") +
    theme_cowplot() +
    theme(legend.position = "top")

# ks test
tpm.melt.all[x.status == "XaXa",mean(value), by=c("gene","chr")][!chr %in% c("chrY",NA), ks.test(V1[chr == "chrX"],V1[chr != "chrX"])]

ggsave2("plots/autosome_ratio_chrX_ecdf.pdf",width = 4,height = 4,
  p.x.ecfd
)

## test expression cut-offs
tpm.auto.ratio.test <- tpm.melt.all[, setNames(lapply(c(0,1,3,10,30,100,300), function(x) median(value[chrx == T & value > x]/median(value[chrx == F & value > x],na.rm = T), na.rm = T) ),paste0(">",c(0,1,3,10,30,100,300))), by=c("sample", "lineage","x.status.simple")]
tpm.auto.ratio.test.melt <- melt(tpm.auto.ratio.test, id.vars= c("sample","lineage","x.status.simple"))

p.auto.ratio.test <- 
  ggplot(tpm.auto.ratio.test.melt[x.status.simple == "XaXa"], aes(y=value, x=variable, col="EPI XaXa")) +
  geom_boxplot() +
  geom_hline(yintercept=1, col=pal.t10[8], lty=2) +
  labs(y="chrX:Autosome ratio", x="Threshold (TPM)") +
  scale_color_brewer(palette="Set1", direction = -1) +
  coord_cartesian(ylim=c(0,1.5)) +
  theme_cowplot() +
  theme(legend.position = "top")

ggsave2("plots/autosome_ratio_cutoff_test.pdf",width = 4,height = 4,
  p.auto.ratio.test
)

suppressMessages(source("R/load_data_cellrep.R"))

## chr:Autosome expression ratios

# load escapees
x.escapees <- fread("x_escapees.tsv")

tpm.melt.all[value > 1 & !gene %in% x.escapees$gene, ratio1 := value/median(value[chrx == F],na.rm = T), by=c("sample", "lineage","x.status.simple")]
tpm.melt.all[value > 10 & !gene %in% x.escapees$gene, ratio10 := value/median(value[chrx == F],na.rm = T), by=c("sample", "lineage","x.status.simple")]
tpm.auto.ratio <- tpm.melt.all[,lapply(.SD,median,na.rm = T),by=c("sample","lineage","sex","chr","x.status.simple"), .SDcols=c("ratio1", "ratio10")] %>% .[,Chromosome := ifelse(chr == "chrX", "chrX",ifelse(chr == "chrY", "chrY","Autosome"))]
# set female chrY NA values to 0 to avoid factor dropping
tpm.auto.ratio[Chromosome == "chrY" & sex == "F" & is.na(ratio1), c("ratio1","ratio10") := 0]

# chr:Autosome bootstrap ratios
# cutoff > 1 tpm
tpm.boot.melt <- as.data.table(melt(ratio.boot(tpm.all,refseq.gene$chrom)))
# cutoff >10 tpm
tpm.boot.melt$value10 <- melt(ratio.boot(tpm.all,refseq.gene$chrom,cutoff = 10))$value
# add metadata
tpm.boot.melt %<>% .[meta,,on="Var2 == LibraryName"]
tpm.boot.melt[Var1 == "chrY" & Sex == "F" & is.na(value), value := 0]
colnames(tpm.boot.melt)[-2] <- c("Chromosome",tolower(colnames(tpm.boot.melt)[-1:-2]))

# plot
library(ggplot2)
library(cowplot)
# create plot skeletons
p.auto.ratio.skel.epi <- 
  ggplot(tpm.auto.ratio[!is.na(Chromosome) & lineage == "EPI"],aes(x=paste(sex, x.status.simple), col=Chromosome)) + 
    geom_hline(yintercept = c(1,0.5), col=pal.t10[8], lty=2) + 
    facet_wrap(~lineage) +
    labs(y="chr:Autosome ratio") + 
    scale_color_brewer(drop=F,palette="Set1",direction = -1) + 
    coord_cartesian(ylim=c(0,1.5)) + 
    cowplot::theme_cowplot() + 
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.auto.ratio.skel.extra <- 
  ggplot(tpm.auto.ratio[!is.na(Chromosome) & lineage != "EPI"],aes(x=sex, col=Chromosome)) + 
  geom_hline(yintercept = c(1,0.5), col=pal.t10[8], lty=2) + 
  facet_wrap(~lineage) +
  labs(y="chr:Autosome ratio") + 
  scale_color_brewer(drop=F,palette="Set1", direction = -1) + 
  coord_cartesian(ylim=c(0,1.5)) + 
  cowplot::theme_cowplot() + 
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

# tpm 1 cutoff
ggsave("plots/autosome_ratio_cutoff_1.pdf",width = 12,height = 4,
plot_grid(nrow = 1,
  p.auto.ratio.skel.epi + geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1, aes(y=ratio1)),
  p.auto.ratio.skel.epi + stat_summary(fun.data = "median_cl_boot",position = position_dodge(0.75),shape=18,stroke=NA,data = tpm.boot.melt[lineage == "EPI"],aes(y=value)),
  p.auto.ratio.skel.extra + geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1, aes(y=ratio1)),
  p.auto.ratio.skel.extra + stat_summary(fun.data = "median_cl_boot",position = position_dodge(0.75),shape=18,stroke=NA,data = tpm.boot.melt[lineage != "EPI"],aes(y=value))
)
)
# tpm 10 cutoff
ggsave("plots/autosome_ratio_cutoff_10.pdf",width = 12,height = 4,
plot_grid(nrow = 1,
  p.auto.ratio.skel.epi + geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1, aes(y=ratio10)),
  p.auto.ratio.skel.epi + stat_summary(fun.data = "median_cl_boot",position = position_dodge(0.75),shape=18,stroke=NA,data = tpm.boot.melt[lineage == "EPI"],aes(y=value10)),
  p.auto.ratio.skel.extra + geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1, aes(y=ratio10)),
  p.auto.ratio.skel.extra + stat_summary(fun.data = "median_cl_boot",position = position_dodge(0.75),shape=18,stroke=NA,data = tpm.boot.melt[lineage != "EPI"],aes(y=value10))
)
)