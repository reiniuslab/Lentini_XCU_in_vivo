suppressMessages(source("R/load_data_smartseq3.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

tpm.melt.ss3[gene.anno.ss3, pos := start, on="var1 == gene_id"]

tpm.melt.gnavg.ss3 <- tpm.melt.ss3[expressed == T & chrx==T & (day != 0 | sex == "M") & x.status %in% c("XaXa","XaXi","XiXa"),mean(value,na.rm=T), by=c("l1","gene","x.status","pos","sex")]

library(ggplot2)
library(cowplot)
p.xexprs.genome.ss3 <-
  ggplot(tpm.melt.gnavg.ss3[!gene %in% x.escapees$gene], aes(y=log2(V1+1), x=pos/1e6, col=l1, fill=l1)) +
    geom_point(alpha=0.1, stroke=NA) +
    geom_smooth() +
    facet_wrap(~sex+x.status, nrow = 1) +
    coord_cartesian(ylim=c(0,6)) +
    labs(y="log2(TPM+1)", x="Chromosome X (Mbp)") +
    scale_color_brewer(palette="Set1", name="XmXp",labels=c("Maternal", "Paternal"), aesthetics = c("col","fill")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_genome.pdf",width = 5,height = 3, p.xexprs.genome.ss3)