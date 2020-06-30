suppressMessages(source("R/load_data_science_ss1.R"))
suppressMessages(source("R/load_data_science.R"))
suppressMessages(source("R/load_data_cellrep.R"))

##

# get X expression fractions
meta.early$c57.x.frac <- get.x.frac.ref(counts.c57.early, counts.cast.early, refseq.gene.early$chrom, prob=1)
meta.blast$c57.x.frac <- get.x.frac.ref(counts.c57.blast, counts.cast.blast, refseq.gene.blast$chrom, prob=1)

tpm.melt.early[meta.early, c57.x.frac := c57.x.frac, on = "sample == LibraryName"]
tpm.melt.blast[meta.blast, c57.x.frac := c57.x.frac, on = "sample == LibraryName"]

# calculate cell averages
tpm.melt.early.avg <- tpm.melt.early[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx","maternal","c57.x.frac")]
tpm.melt.early.avg[, maternal.x.frac := ifelse(maternal == "C57",c57.x.frac,1-c57.x.frac)]

tpm.melt.blast.avg <- tpm.melt.blast[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx","maternal","c57.x.frac")]
tpm.melt.blast.avg[, maternal.x.frac := ifelse(maternal == "C57",c57.x.frac,1-c57.x.frac)]

tpm.melt.merge.avg <- rbindlist(list(tpm.melt.early.avg,tpm.melt.blast.avg))

# plot
library(ggplot2)
library(cowplot)
p.exprs.m.frac.preimpl <- 
  ggplot(tpm.melt.merge.avg[is.na(sex) | sex == "F"],aes(y = V1, x = maternal.x.frac, col = paste(chrx,l1==tolower(maternal)), fill = paste(chrx,l1==tolower(maternal)) )) +
    geom_vline(xintercept = 0.5, col=pal.t10[8], lty=2) +
    geom_point(stroke=NA, alpha=0.1) +
    geom_smooth(method="lm") +
    labs(x = "Maternal X-fraction", y = "Expression (TPM)") +
    coord_cartesian(xlim=c(0,1), ylim=c(0,100)) +
    expand_limits(y=0) +
    facet_grid(lineage ~ chrx) +
    scale_fill_brewer(palette="Paired", name = NULL, labels =c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal")) +
    scale_color_brewer(palette = "Paired", name = NULL, labels =c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/m_fraction_expression_preimpl.pdf", width = 4, height = 12,
  p.exprs.m.frac.preimpl
)