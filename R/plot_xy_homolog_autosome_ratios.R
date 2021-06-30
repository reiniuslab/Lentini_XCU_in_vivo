suppressMessages(source("R/load_data_cellrep.R"))

# add homology information
# x-y homologs modified from Table1: doi:10.1016/j.cell.2014.09.052
xy_hom <- fread("xy_homologs.tsv")

tpm.melt.all[,xy_homolog := gene %in% xy_hom[,c(Ygene,Xhomolog)]]
tpm.melt.all[xy_hom, c("class", "group") := list(SequenceClass, Group), on=c("gene == Ygene")]
tpm.melt.all[xy_hom, c("class", "group") := list(SequenceClass, Group), on=c("gene == Xhomolog")]
tpm.melt.all[xy_hom[,paste0(unique(c(Xhomolog,Ygene)),collapse = ":"),by="Group"], group_gene := V1, on=c("group == Group")]

# calculate X:A ratios based on XY homology
tpm.melt.all[value > 1, ratio1 := value/median(value[chrx == F],na.rm = T), by=c("sample", "lineage","x.status.simple")]
tpm.auto.ratio.xy <- tpm.melt.all[,median(ratio1,na.rm = T),by=c("sample","lineage","sex","chr","x.status.simple","xy_homolog")] %>% .[,Chromosome := ifelse(chr == "chrX", "chrX",ifelse(chr == "chrY", "chrY","Autosome"))]
# set chrY NA values to 0 to avoid factor dropping
tpm.auto.ratio.xy[Chromosome == "chrY" & is.na(V1), V1 := 0]

# P-values vs. expected
tpm.auto.ratio.xy[!is.na(Chromosome) & chr == "chrX" & xy_homolog == T, wilcox.test(V1, mu = 0.5)$p.value, by=c("lineage", "sex")]

## PLOT
library(ggplot2)
library(cowplot)
# check expressed chrY genes
p.det.y <- 
  ggplot(tpm.melt.all[chr == "chrY", mean(value > 0), by = c("gene","sex")], aes(y= V1, x = reorder(gene,V1), fill = sex)) +
  geom_col(position = position_dodge(0.9)) +
  labs(y="fraction detected chrY genes", x = NULL) +
  coord_flip() +
  scale_fill_brewer(palette="Set1") +
  theme_cowplot()

# Autosomal ratios
p.auto.ratio.xy.epi <-
  ggplot(tpm.auto.ratio.xy[!is.na(Chromosome) & lineage == "EPI" & chr %in% c("chrY","chrX")],aes(y = V1, x=paste(sex, x.status.simple), col=paste(chr,xy_homolog) )) + 
    geom_hline(yintercept = c(1,0.5), col=pal.t10[8], lty=2) + 
    geom_boxplot(outlier.alpha = 0.1,outlier.stroke = NA) +
    facet_wrap(~lineage) +
    labs(y="chr:Autosome ratio") + 
    scale_color_brewer(name=NULL,palette="Set1",direction = -1,labels=c("chrX:nonHom","chrx:hom","chrY:nonHom","chrY:hom")) + 
    coord_cartesian(ylim=c(0,1.5)) + 
    cowplot::theme_cowplot() + 
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.auto.ratio.xy.extra <-
  ggplot(tpm.auto.ratio.xy[!is.na(Chromosome) & lineage != "EPI" & chr %in% c("chrY","chrX")],aes(y = V1, x=sex, col=paste(chr,xy_homolog) )) + 
    geom_hline(yintercept = c(1,0.5), col=pal.t10[8], lty=2) + 
    geom_boxplot(outlier.alpha = 0.1,outlier.stroke = NA) +
    facet_wrap(~lineage) +
    labs(y="chr:Autosome ratio") + 
    scale_color_brewer(name=NULL,palette="Set1",direction = -1,labels=c("chrX:nonHom","chrx:hom","chrY:nonHom","chrY:hom")) + 
    coord_cartesian(ylim=c(0,1.5)) + 
    cowplot::theme_cowplot() + 
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave("plots/xy_homologs_autosome_ratios.pdf",height = 4,width = 6,
  plot_grid(nrow=1,
    p.auto.ratio.xy.epi,
    p.auto.ratio.xy.extra
  )
)