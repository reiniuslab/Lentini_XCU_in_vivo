suppressMessages(source("R/load_data_cellrep.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

# non-allele resolution
tpm.melt.all[expressed_lineage == T & !gene %in% x.escapees$gene, mf_ratio := value / mean(value[sex == "M"], na.rm = T), by=c("gene","embryonicday","lineage")]
tpm.melt.all.gene.avg <- tpm.melt.all[, mean(mf_ratio,na.rm = T), by=c("gene","x.status.maternal","lineage","sex","chrx")]

# allele resolution
tpm.melt[, l1_x := ifelse(l1 == tolower(maternal),substr(x.status.maternal,1,2),substr(x.status.maternal,3,4))]
tpm.melt[, activex := ifelse(l1_x == "Xa",T,F)]
tpm.melt[expressed_lineage == T & !gene %in% x.escapees$gene, mf_ratio := value / mean(value[sex == "M"], na.rm = T), by=c("gene","embryonicday","lineage","activex")]
## fix XaXa alleles
tpm.melt[,allele := 1]
tpm.melt[x.status.maternal == "XaXa" & chrx == T, allele := ifelse(l1 == tolower(maternal),1,2) ]
tpm.melt.gene.avg <- tpm.melt[,mean(mf_ratio,na.rm = T), by=c("gene","x.status.maternal","lineage","sex","l1","chrx","activex","allele")]

# plot
library(ggplot2)
library(cowplot)
p.mf.ratio.all.epi <-  
  ggplot(tpm.melt.all.gene.avg[!is.na(chrx) & lineage == "EPI" & sex == "F"], aes(y = log2(V1), x = factor(x.status.maternal,levels = c("XiXa","XsXa","XaXa","XaXs","XaXi")), col=chrx )) +
    stat_summary(fun.data="median_cl_boot", shape=18, stroke=NA) +
    geom_hline(yintercept = 0, col=pal.t10[8], lty=2) +
    ylab("log2 F:M ratio") +
    coord_cartesian(ylim=c(-1,1)) +
    facet_grid(~lineage, scales = "free_x",space = "free_x") +
    scale_colour_brewer(palette="Set1") +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.title.x = element_blank())

p.mf.ratio.all.extra <- 
  ggplot(tpm.melt.all.gene.avg[!is.na(chrx) & lineage != "EPI" & sex == "F"], aes(y = log2(V1), x = sex, col=chrx )) +
    stat_summary(fun.data="median_cl_boot", shape=18, stroke=NA) +
    geom_hline(yintercept = 0, col=pal.t10[8], lty=2) +
    coord_cartesian(ylim=c(-1,1)) +
    facet_grid(~lineage, scales = "free_x",space = "free_x") +
    scale_colour_brewer(palette="Set1") +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y =element_blank(), axis.text.y = element_blank())

p.mf.ratio.allele.epi <- 
  ggplot(tpm.melt.gene.avg[!is.na(chrx) & lineage == "EPI" & sex == "F"], aes(y = log2(V1), x = factor(x.status.maternal,levels = c("XiXa","XsXa","XaXa","XaXs","XaXi")), x2= factor(allele), col=chrx )) +
    stat_summary(fun.data="median_cl_boot", shape=20, stroke=NA) +
    geom_hline(yintercept = 0, col=pal.t10[8], lty=2) +
    ylab("log2 F:M ratio") +
    coord_cartesian(ylim=c(-1,1)) +
    facet_grid(~lineage, scales = "free_x",space = "free_x") +
    scale_colour_brewer(palette="Set1") +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.title.x = element_blank())

p.mf.ratio.allele.extra <- 
  ggplot(tpm.melt.gene.avg[!is.na(chrx) & lineage != "EPI" & sex == "F"], aes(y = log2(V1), x = sex, col=chrx )) +
    stat_summary(fun.data="median_cl_boot", shape=20, stroke=NA) +
    geom_hline(yintercept = 0, col=pal.t10[8], lty=2) +
    coord_cartesian(ylim=c(-2,2)) +
    facet_grid(~lineage, scales = "free_x",space = "free_x") +
    scale_colour_brewer(palette="Set1") +
    theme_cowplot() +
    theme(legend.position = "top",strip.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y =element_blank(), axis.text.y = element_blank())

ggsave2("plots/mf_ratio_allele_vs_non.pdf",width = 4,height = 5,
  plot_grid(nrow = 2,rel_widths = c(1,0.45),
    p.mf.ratio.all.epi, p.mf.ratio.all.extra,
    p.mf.ratio.allele.epi, p.mf.ratio.allele.extra
  )
)