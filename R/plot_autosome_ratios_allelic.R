suppressMessages(source("R/load_data_cellrep.R"))

## Allelic chr:Autosome expression ratios

# load escapees
x.escapees <- fread("x_escapees.tsv")

tpm.melt[value > 1 & !gene %in% x.escapees$gene, ratio := value/median(value[chrx == F],na.rm = T), by=c("sample", "lineage","x.status", "l1")]
tpm.auto.ratio.allele <- tpm.melt[,median(ratio,na.rm = T),by=c("sample","lineage","sex","chrx","x.status", "l1")]

# plot
library(ggplot2)
library(cowplot)
# EPI
p.auto.ratio.allele.epi <-
ggplot(tpm.auto.ratio.allele[lineage == "EPI" & chrx == T], aes(y=V1, x=factor(x.status,levels=c("XiXa","XsXa","XaXa","XaXs","XaXi")), col=l1) ) +
  geom_hline(yintercept = c(1,0.5), col=pal.t10[8], lty=2) + 
  geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
  facet_grid(~lineage+sex, scales = "free_x", space="free_x") +
  labs(y="Allelic X:Autosome ratio") + 
  scale_color_brewer(drop=F,palette="Set1",direction = -1, name="XC57XCAST") + 
  coord_cartesian(ylim=c(0,1.5)) + 
  cowplot::theme_cowplot() + 
  theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

# save
ggsave2("plots/autosome_ratio_allelic.pdf",width = 4,height = 4,
  p.auto.ratio.allele.epi
)