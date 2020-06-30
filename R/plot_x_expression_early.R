suppressMessages(source("R/load_data_science_ss1.R"))

# average
tpm.melt.early.avg <- tpm.melt.early[expressed_lineage == T & chr != "chrY", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","lineage","sex","l1","chrx","maternal")]

## plot
library(ggplot2)
library(cowplot)
p.x.expr.early <- 
  ggplot(tpm.melt.early.avg, aes(y = V1, x= lineage, x2= l1, col=paste(chrx, tolower(maternal) == l1))) +
    geom_boxplot(outlier.alpha = 0.1, outlier.stroke = NA) +
    coord_cartesian(ylim=c(0,100)) +
    ylab("Expression (TPM)") +
    facet_wrap(~sex) +
    scale_color_brewer(name = "C57xCAST:",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/x_expression_boxplot_early.pdf",height = 4, width = 10,
  p.x.expr.early
)