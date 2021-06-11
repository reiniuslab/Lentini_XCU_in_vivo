suppressMessages(source("R/load_data_cellrep.R"))

# average
tpm.melt.avg <- tpm.melt[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample","lineage","sex","l1","embryonicday","x.status.maternal","chrx", "maternal")]
tpm.melt.all.avg <- tpm.melt.all[expressed_lineage == T & chrx == T, mean(value, na.rm=T), by=c("sample", "lineage", "sex", "embryonicday", "x.status.simple", "c57.x.frac")]

## plot
library(ggplot2)
library(cowplot)
# by x status
p.x.expr.tpm.allele.epi <- 
  ggplot(tpm.melt.avg[lineage == "EPI"],aes(y=V1,x= factor(paste(sex,x.status.maternal), levels = paste(rep(c("F","M"),each=5),c("XaXi","XaXs","XaXa","XsXa","XiXa")) ),col=paste(chrx, tolower(maternal) == l1))) +
  geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
  geom_vline(xintercept = seq(1.5,5.5), col=pal.t10[8], lty=2) +
  expand_limits(y=0) +
  ylab("Expression (TPM)") +
  coord_cartesian(ylim=c(0,60)) +
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
  coord_cartesian(ylim=c(0,60)) +
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

# by e day
p.x.expr.tpm.eday <- 
  ggplot(tpm.melt.avg[lineage == "EPI" & chrx == T],aes(y=V1,x=factor(embryonicday),col=tolower(maternal) == l1) ) +
  stat_summary(fun.data="median_cl_boot", stroke=NA) +
  expand_limits(y=0) +
  labs(y="Expression (TPM)", x="Embryonic day") +
  coord_cartesian(ylim=c(0,60)) +
  facet_grid(~sex+x.status.maternal) +
  scale_color_brewer(name = "ChrX",labels = c("Paternal","Maternal"), palette = "Paired") +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/x_expression_eday.pdf",height = 3, width = 6, p.x.expr.tpm.eday)

# total expression
p.x.expr.tpm.total.epi <-
ggplot(tpm.melt.all.avg[lineage == "EPI"], aes(x=x.status.simple, y=V1)) + 
  geom_boxplot(outlier.stroke = NA) + 
  facet_grid(~lineage+sex, scales="free_x", space="free_x") +
  labs(x=NULL, y="Expression (TPM)") +
  expand_limits(y=0) +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/x_expression_total.pdf",height = 3, width = 3,
  p.x.expr.tpm.total.epi
)
