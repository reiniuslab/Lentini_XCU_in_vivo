suppressMessages(source("R/load_data_cellrep.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

# add gene position
tpm.melt[refseq.gene, pos := txStart, on = "gene==name2"]

# tpm.melt.median <- tpm.melt[chrx == T & gene != x.esc.gn & expressed_lineage == T, mean(value, na.rm = T), by= c("lineage", "sex", "x.status","l1","gene")] %>% .[,log(median(V1, na.rm = T)+1), by=c("lineage","sex","x.status","l1")]

library(ggplot2)
library(cowplot)
p.xexprs.genome.epi <-
  ggplot(tpm.melt[lineage == "EPI" & chrx == T & !gene %in% x.escapees$gene & expressed_lineage == T], aes(y = log2(value+1), x = pos/1e6, col = l1 == tolower(maternal), fill= l1 == tolower(maternal) )) +
    stat_summary(fun.y = "mean", geom="point", alpha=0.1, stroke=NA, show.legend = F) +
    geom_smooth() +
    coord_cartesian(ylim = c(0,7)) +
    facet_wrap(~lineage + sex + x.status.maternal, scales = "free_x",nrow = 1) +
    labs(x="Chromosome X (Mbp)", y="log2(TPM + 1)") +
    scale_color_brewer(name = "Allele (XMXP)",labels = c("Paternal","Maternal"), palette="Set1") +
    scale_fill_brewer(name = "Allele (XMXP)",labels = c("Paternal","Maternal"), palette="Set1") +
    guides(fill=F) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

p.xexprs.genome.extra <-
  ggplot(tpm.melt[lineage != "EPI" & chrx == T & !gene %in% x.escapees$gene & expressed_lineage == T], aes(y = log2(value+1), x = pos/1e6, col = l1 == tolower(maternal), fill= l1 == tolower(maternal) )) +
    stat_summary(fun.y = "mean", geom="point", alpha=0.1, stroke=NA, show.legend = F) +
    geom_smooth() +
    coord_cartesian(ylim = c(0,7)) +
    facet_wrap(~lineage + sex + 1, scales = "free_x",nrow = 1) +
    labs(x="Chromosome X (Mbp)", y="log2(TPM + 1)") +
    scale_color_brewer(name = NULL,labels = c("Paternal","Maternal"), palette="Set1") +
    scale_fill_brewer(name = NULL,labels = c("Paternal","Maternal"), palette="Set1") +
    guides(fill=F) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/x_expression_genome.pdf", width = 12, height = 3,
  plot_grid(nrow = 1, rel_widths = c(1,0.66),
    p.xexprs.genome.epi, 
    p.xexprs.genome.extra
  )
)