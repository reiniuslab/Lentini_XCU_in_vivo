suppressMessages(source("R/load_data_smartseq3.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

# add gene positions
si <- SeqinfoForBSGenome("mm10")

tpm.melt.ss3[gene.anno.ss3, pos := start, on="var1 == gene_id"]
tpm.melt.ss3[, pos_max := seqlengths(si)[chr]]

tpm.melt.gnavg.ss3 <- tpm.melt.ss3[expressed == T & x.status %in% c("XaXa","XaXi","XiXa"),mean(value,na.rm=T), by=c("l1","gene","chr","x.status","pos","pos_max","sex")]
tpm.melt.gnavg.ss3[, chr := factor(chr, levels=paste0("chr", c(1:19,"X")))]

library(ggplot2)
library(cowplot)
# chr4 + chrX
p.xexprs.genome.ss3 <-
  ggplot(tpm.melt.gnavg.ss3[!gene %in% x.escapees$gene & chr %in% c("chr4","chrX")], aes(y=log2(V1+1), x=pos/pos_max, col=paste(chr,l1=="c57"), fill=paste(chr,l1=="c57") )) +
    geom_smooth(method="loess") +
    facet_grid(~x.status) +
    coord_cartesian(ylim=c(0,6)) +
    labs(y="log2(TPM+1)", x="Chromosome (%)") +
    scale_color_brewer(palette="Paired", name=NULL, labels=c("Chr4:Paternal", "Chr4:Maternal", "ChrX:Paternal", "ChrX:Maternal"), aesthetics = c("col","fill")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_genome_new.pdf",width = 5,height = 3, p.xexprs.genome.ss3)

# all chr
## by chromosome position
p.xexprs.genome.full.ss3 <-
  ggplot(tpm.melt.gnavg.ss3[!gene %in% x.escapees$gene], aes(y=log2(V1+1), x=pos/1e6, col=l1=="c57", fill=l1=="c57" )) +
  geom_smooth(method="loess") +
  facet_grid(x.status~chr, scales="free_x", space="free_x") +
  coord_cartesian(ylim=c(0,6)) +
  labs(y="log2(TPM+1)", x="Chromosome (Mbp)") +
  scale_color_brewer(palette="Paired", name=NULL, labels=c("Paternal", "Maternal"), aesthetics = c("col","fill")) +
  theme_cowplot() +
  theme(legend.position = "top", strip.background = element_blank())

## autosomal average
p.xexprs.genome.full.avg.ss3 <-
  ggplot(tpm.melt.gnavg.ss3[chr != "chrX"], aes(y=log2(V1+1), x=pos/pos_max, col=l1=="c57", fill=l1=="c57" )) +
  geom_smooth(method="loess") +
  facet_grid(x.status~"global") +
  coord_cartesian(ylim=c(0,6)) +
  labs(y="log2(TPM+1)", x="Chromosome (%)") +
  scale_color_brewer(palette="Paired", name=NULL, labels=c("Paternal", "Maternal"), aesthetics = c("col","fill")) +
  theme_cowplot() +
  theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_genome_allchr.pdf",width = 16,height = 6, 
  plot_grid(nrow = 1, rel_widths = c(9,1), p.xexprs.genome.full.ss3, p.xexprs.genome.full.avg.ss3)
)

# Plot ideogram
library(Gviz)
ideo <- IdeogramTrack(genome = "mm10")
si <- keepStandardChromosomes(SeqinfoForBSGenome(genome = "mm10"))

pdf("plots/ss3_x_expression_genome_ideogram.pdf")
  lapply(seq_along(si),function(i) plotTracks(ideo,chromosome=seqnames(si)[i], main=seqnames(si)[i],from=1,to=seqlengths(si)[i]) )
dev.off()