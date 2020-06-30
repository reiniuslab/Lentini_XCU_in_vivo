suppressMessages(source("R/load_data_cellrep.R"))

# tpm per allele
tpm.melt.avg <- tpm.melt[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample","lineage","sex","l1","embryonicday","x.status.maternal","chrx", "maternal")]

## plot
library(ggplot2)
library(cowplot)
# by embryonic day
p.x.expr.tpm.eday <-
  ggplot(tpm.melt.avg,aes(y=V1,x=factor(paste0("E",embryonicday)),col=paste(chrx, tolower(maternal) == l1))) +
  geom_jitter(stroke = NA, width=0.2) +
  #geom_vline(xintercept = seq(1.5,5.5), col=pal.t10[8], lty=2) +
  expand_limits(y=0) +
  ylab("Expression (TPM)") +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(~sex + lineage, scales = "free_x",space = "free_x") +
  scale_color_brewer(name = "C57xCAST",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/x_expression_scatter_eday.pdf",height = 4, width = 6,
  p.x.expr.tpm.eday
)
# by x status
p.x.expr.tpm.allele.epi <- 
  ggplot(tpm.melt.avg[lineage == "EPI"],aes(y=V1,x= factor(paste(sex,x.status.maternal), levels = paste(rep(c("F","M"),each=5),c("XiXa","XsXa","XaXa","XaXs","XaXi")) ),col=paste(chrx, tolower(maternal) == l1))) +
  geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
  geom_vline(xintercept = seq(1.5,5.5), col=pal.t10[8], lty=2) +
  expand_limits(y=0) +
  ylab("Expression (TPM)") +
  coord_cartesian(ylim=c(0,100)) +
  facet_wrap(~lineage) +
  scale_color_brewer(name = "C57xCAST",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.x.expr.tpm.allele.extra <- 
  ggplot(tpm.melt.avg[lineage != "EPI"],aes(y=V1,x=sex,col=paste(chrx, tolower(maternal) == l1))) +
  geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
  geom_vline(xintercept = 1.5, col=pal.t10[8], lty=2) +
  expand_limits(y=0) +
  ylab("Expression (TPM)") +
  coord_cartesian(ylim=c(0,150)) +
  facet_wrap(~lineage) +
  scale_color_brewer(name = "C57xCAST",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/x_expression_boxplot.pdf",height = 4, width = 10,
  plot_grid(nrow = 1, rel_widths = c(1,0.66),
    p.x.expr.tpm.allele.epi,
    p.x.expr.tpm.allele.extra
  )
)