suppressMessages(source("R/load_data_smartseq3.R"))

# load subsampling data
fls <- c(list.files("smartseq3/results/zUMIs_output/expression",pattern = "mESC_day0.umicount.exon.downsampled",full.names = T),"smartseq3/results/zUMIs_output/expression/mESC_day0.umicount.exon.all.loom")

relumi.subsample.ls <- lapply(fls,function(x){
  dat <- loom2matrix(x)
  dat.rel <- apply(dat,2,umi2rel)
  dat.rel.melt <- as.data.table(melt(dat.rel))
  dat.rel.melt %<>% add.info.ss3(coldata = meta.ss3[day==0])
  dat.rel.melt %<>% .[exclude == F]
  dat.rel.melt[,expressed := mean(value,na.rm = T)>0, by = c("gene")]
  return(dat.rel.melt)
})
names(relumi.subsample.ls) <- sapply(strsplit(basename(fls),"\\."),"[[", 4)

relumi.melt.subsample <- rbindlist(relumi.subsample.ls,idcol=T)
relumi.melt.subsample.avg <- relumi.melt.subsample[expressed == T & chr != "chrY", mean(value,na.rm=T),by=c("sample_id",".id","sex","day","chrx","x.status")]

library(ggplot2)
library(cowplot)
p.ss3.subsample.umi <-
  ggplot(relumi.melt.subsample.avg[x.status %in% c("XaXa","XaXi","XiXa")], aes(y=V1, x=factor(gsub("downsampled_","",.id), levels=c("50000","100000","250000","500000","750000","1000000","all")), col=chrx)) +
    stat_summary(fun.y="median", geom="line", show.legend = F, aes(group=chrx)) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    facet_wrap(~sex+x.status, nrow=1) +
    labs(y="Relative UMIs", x="Subsampled reads (x10^6)") +
    coord_cartesian(ylim=c(0,1)) +
    scale_x_discrete(labels=c("0.05","0.10","0.25","0.50","0.75","1.00","all")) +
    scale_colour_brewer(palette="Set1", name=NULL, labels=c("Autosomes","ChrX")) +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.text.x = element_text(angle=90,vjust=0.5, hjust=1))

p.ss3.subsample.ratio <-
  ggplot(relumi.melt.subsample.avg[,V1[chrx == T] / median(V1[chrx == F], na.rm = T), by=c("sample_id", ".id", "sex", "day", "x.status")][x.status %in% c("XaXa","XaXi","XiXa")], aes(y=V1, x=factor(gsub("downsampled_","",.id), levels=c("50000","100000","250000","500000","750000","1000000","all")),col="Ratio" )) +
    stat_summary(fun.y="median", geom="line", show.legend = F, aes(group = paste(sex,x.status))) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    facet_wrap(~sex+x.status, nrow=1) +
    labs(y="Ratio", x="Subsampled reads (x10^6)") +
    coord_cartesian(ylim=c(0,2)) +
    scale_x_discrete(labels=c("0.05","0.10","0.25","0.50","0.75","1.00","all")) +
    scale_colour_manual(values = "black", name=NULL) +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.text.x = element_text(angle=90,vjust=0.5, hjust=1))

ggsave2("plots/ss3_subsample_x_expression.pdf",width = 4,height = 6,
  plot_grid(ncol = 1,align = "v",
    p.ss3.subsample.umi,
    p.ss3.subsample.ratio
  )
)