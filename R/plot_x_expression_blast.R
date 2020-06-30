suppressMessages(source("R/load_data_science.R"))

# average
tpm.melt.blast.avg <- tpm.melt.blast[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx", "maternal")]
tpm.melt.blast.avg[,lineage := factor(lineage,levels=c("earlyblast","midblast","lateblast"))]

## plot
library(ggplot2)
library(cowplot)
p.x.expr.blast <- 
  ggplot(tpm.melt.blast.avg,aes(y=V1,x=lineage,x2=factor(l1,levels=c("cast","c57")),col=paste(chrx, tolower(maternal) == l1))) +
    geom_boxplot(outlier.alpha = 0.1, outlier.stroke = NA) +
    expand_limits(y=0) +
    ylab("Expression (TPM)") +
    coord_cartesian(ylim=c(0,100)) +
    facet_wrap(~sex) +
    scale_color_brewer(name = "C57xCAST",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/x_expression_boxplot_blast.pdf",height = 4, width = 6,
  p.x.expr.blast
)

## by subcluster
tpm.melt.blast.avg %<>% .[fread("data/blast_subcluster.tsv"), on="sample == V1"]
tpm.melt.blast.avg[, label := cluster] %>% .[!is.na(cluster2), label := cluster2]

p.x.expr.blast.sub <- 
  ggplot(tpm.melt.blast.avg[label %in% c(2,5,6) & lineage == "lateblast"],aes(y=V1,x=paste(label,sex),x2=factor(l1,levels=c("cast","c57")),col=paste(chrx, tolower(maternal) == l1))) +
  geom_boxplot(outlier.alpha = 0.1, outlier.stroke = NA) +
  expand_limits(y=0) +
  ylab("Expression (TPM)") +
  coord_flip(ylim=c(0,150)) +
  facet_wrap(~lineage) +
  scale_color_brewer(name = "C57xCAST:",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top", axis.title.y = element_blank(), strip.background = element_blank())

ggsave2("plots/x_expression_boxplot_blast_subcluster.pdf",height = 6, width = 3,
  p.x.expr.blast.sub
)