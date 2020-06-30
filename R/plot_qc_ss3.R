suppressMessages(source("R/load_data_smartseq3.R"))

library(ggplot2)
library(cowplot)

# read depth
stats.rpc.ss3 <- rbindlist(idcol=T,
  list(
    "day0" = fread("smartseq3/results/zUMIs_output/stats/mESC_day0.readspercell.txt"),
    "cl5diff" = fread("smartseq3/results/zUMIs_output/stats/mESC_diff_clone5.readspercell.txt") 
  )
)
stats.rpc.ss3[meta.ss3, exclude := exclude, on="RG == sample_bc"]

# umi counts
stats.umi.ss3 <- rbindlist(idcol=T,
  list(
    "day0" = fread("smartseq3/results/zUMIs_output/stats/mESC_day0.UMIcounts.txt"),
    "cl5diff" = fread("smartseq3/results/zUMIs_output/stats/mESC_diff_clone5.UMIcounts.txt") 
  )
)
stats.umi.ss3[meta.ss3, exclude := exclude, on="SampleID == sample_bc"]

# detected genes
gn.stat <- tpm.melt.all.ss3[,list(all=sum(value > 0)),by="sample_id"][tpm.melt.ss3[,sum(value > 0,na.rm = T),by=c("sample_id","l1")][,max(V1),by="sample_id"], allele := V1, on = "sample_id == sample_id"]
gn.stat[,ratio := allele/all]

ggsave2("plots/ss3_qc_stats.pdf",width = 4,height = 3,
  plot_grid(ncol=2,
    ggplot(stats.rpc.ss3[exclude == F & type == "Exon"], aes(x=N/1e6)) + geom_histogram() + labs(y="Count", x="Read depth (x1,000,000)") + expand_limits(x=0) + theme_cowplot(),
    ggplot(gn.stat, aes(x=all/1e3)) + labs(y="Count",x="Genes (x1,000)") + geom_histogram() + expand_limits(x=0) + theme_cowplot(),
    ggplot(stats.umi.ss3[exclude == F & type == "Exon"], aes(x=Count/1e3)) + geom_histogram() + labs(y="Count", x="UMI counts (x1,000)") + expand_limits(x=0) + theme_cowplot(),
    ggplot(gn.stat, aes(x=allele/1e3)) + labs(y="Count",x="Allele-level (x1,000)") + geom_histogram() + expand_limits(x=0) + theme_cowplot(),
    qplot(stats.umi.ss3[exclude == F & type == "Exon",Count]/stats.rpc.ss3[exclude == F & type == "Exon"][order(-.id,RG),N]*100) + labs(y="Count",x="%UMI reads") + theme_cowplot(),
    ggplot(gn.stat, aes(x=ratio*100)) + labs(y="Count",x="%allelic") + geom_histogram() + expand_limits(x=0) + theme_cowplot()
  )
)

library(ggbeeswarm)
p.markers.ss3 <- 
  ggplot(tpm.melt.all.ss3[gene %in% c("Fgf5","Sox2","Pou5f1","Nanog")],aes(y=value,x=factor(day) )) +
    geom_quasirandom(stroke=NA, alpha=0.1) +
    geom_violin(scale="width", fill=NA) +
    facet_grid(gene~factor(sex,levels=c("M","F")), scales="free",space="free_x") +
    labs(y="Expression (TPM)", x="Day") +
    theme_cowplot() +
    theme(strip.background = element_blank())

ggsave2("plots/ss3_marker_gene.pdf",width = 3,height = 3, p.markers.ss3)
#

p.xist.ss3 <-
  ggplot(tpm.melt.all.ss3[gene == "Xist"],aes(y=value,x=factor(day) )) +
    geom_quasirandom(stroke=NA, alpha=0.1) +
    geom_violin(scale="width", fill=NA) +
    facet_grid(gene~factor(sex,levels=c("M","F")), scales="free",space="free_x") +
    labs(y="Expression (TPM)", x="Day") +
    theme_cowplot() +
    theme(strip.background = element_blank())

ggsave2("plots/ss3_xist_expression.pdf",width = 3,height = 2, p.xist.ss3)