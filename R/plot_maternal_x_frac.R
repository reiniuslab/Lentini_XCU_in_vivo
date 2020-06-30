suppressMessages(source("R/load_data_cellrep.R"))
suppressMessages(source("R/load_data_science.R"))
suppressMessages(source("R/load_data_science_ss1.R"))

# determine C57 fraction
meta.blast$c57.x.frac <- get.x.frac.ref(counts.c57.blast, counts.cast.blast, refseq.gene.blast$chrom)
meta.early$c57.x.frac <- get.x.frac.ref(counts.c57.early, counts.cast.early, refseq.gene.early$chrom)

# maternal fraction
meta[,maternal.x.frac := ifelse(Maternal == "CAST", 1 - c57.x.frac, c57.x.frac)]
meta.blast[,maternal.x.frac := ifelse(Maternal == "CAST", 1 - c57.x.frac, c57.x.frac)]
meta.early[,maternal.x.frac := ifelse(Maternal == "CAST", 1 - c57.x.frac, c57.x.frac)]

# blastocyst subclusters
meta.blast %<>% .[fread("data/blast_subcluster.tsv"), on="LibraryName == V1"]
meta.blast[, label := cluster] %>% .[!is.na(cluster2), label := cluster2]
meta.blast$label[!meta.blast$label %in% c(2,5,6)] <- "ND"

library(ggplot2)
library(ggbeeswarm)
library(cowplot)
# early development
p.mat.x.f.early <- ggplot(meta.early[Sex %in% c("F",NA)],aes(x=Lineage, y=maternal.x.frac, col="black")) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(dodge.width = 0.9, stroke=NA, show.legend = T) + coord_cartesian(ylim=c(0,1)) + scale_x_discrete(limits=levels(meta.early$Lineage)) + scale_color_manual(values="black") + labs(y="Maternal chrX fraction", x="Early development") + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", strip.background = element_blank())
p.mat.x.m.early <- ggplot(meta.early[Sex == "M"],aes(x="Early", y=maternal.x.frac, col="black")) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(dodge.width = 0.9, stroke=NA) + coord_cartesian(ylim=c(0,1)) + labs(y="Maternal chrX fraction", x="Early development")+ scale_color_manual(values="black") + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())
# blastocyst
p.mat.x.f.blast <- ggplot(meta.blast[Sex == "F"],aes(x=Lineage, y=maternal.x.frac, col=factor(label))) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(na.rm = F, dodge.width = 0.9, stroke=NA) + geom_vline(xintercept = c(1.5,2.5), col=pal.t10[8], lty=2) + coord_cartesian(ylim=c(0,1)) + scale_x_discrete(limits=c("earlyblast","midblast","lateblast")) + labs(y="Maternal chrX fraction", x="Blastocyst") + scale_color_manual(name="Lineage",values=c(pal.t10[c(2,5,6)],"black")) + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())
p.mat.x.m.blast <- ggplot(meta.blast[Sex == "M"],aes(x="Blast", y=maternal.x.frac, col=factor(label))) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(na.rm = F, dodge.width = 0.9, stroke=NA) + coord_cartesian(ylim=c(0,1)) + labs(y="Maternal chrX fraction", x="Blastocyst") + scale_color_manual(name="Lineage",values=c(pal.t10[c(2,5,6)],"black")) + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())
# post-implantation
p.mat.x.f <- ggplot(meta[Sex == "F"],aes(x=factor(EmbryonicDay), y=maternal.x.frac, col=Lineage)) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(dodge.width = 0.9, stroke=NA) + geom_vline(xintercept = c(1.5,2.5), col=pal.t10[8], lty=2) + coord_cartesian(ylim=c(0,1)) + labs(y="Maternal chrX fraction", x="Post-implantation (embryonic day)") + scale_color_manual(values=pal.t10[c(1,3,4)]) + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())
p.mat.x.m <- ggplot(meta[Sex == "M"],aes(x="Post-implantation", y=maternal.x.frac, col=Lineage)) + geom_violin(scale="width",lty=2, aes(fill=NULL)) + geom_quasirandom(dodge.width = 0.9, stroke=NA) + coord_cartesian(ylim=c(0,1)) + labs(y="Maternal chrX fraction", x="Post-implantation") + scale_color_manual(values=pal.t10[c(1,3,4)]) + facet_wrap(~Sex) + theme_cowplot() + theme(legend.position = "top", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())

# plot
ggsave2("plots/maternal_x_fraction.pdf",width = 8,height = 3,
plot_grid(
  nrow = 1,
  rel_widths = c(0.66,0.33,0.33,0.05,0.15,0.1),
  plotlist = list(
    p.mat.x.f.early,
    p.mat.x.f.blast,
    p.mat.x.f,
    p.mat.x.m.early,
    p.mat.x.m.blast,
    p.mat.x.m
))
)