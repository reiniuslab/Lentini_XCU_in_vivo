suppressMessages(source("R/load_data_smartseq3.R"))

altratio.melt.ss3 <- rbindlist(list(
  melt(fread("smartseq3/results/zUMIs_output/allelic/mESC_day0.fract_CAST_reads.txt"),id.vars="GeneID"),
  melt(fread("smartseq3/results/zUMIs_output/allelic/mESC_diff_clone5.fract_CAST_reads.txt"),id.vars="GeneID")
))

altratio.melt.ss3[gene.anno.ss3, chr := chrom, on = "GeneID == gene_id"]
altratio.melt.ss3[meta.ss3, c("day", "sex", "exclude","x.status") := list(day, sex, exclude, x.status), on="variable == sample_bc"]
altratio.melt.ss3[,growth := ifelse(day == 0, "2i", "Diff")]
altratio.melt.ss3[exclude == F, detected := mean(!is.na(value))>0.1, by="GeneID"]
altratio.melt.ss3[, chr := factor(chr, levels=c(1:19,"X"))] # fix chromosome order
altratio.melt.ss3[gene.anno.ss3, pos := start, on = "GeneID == gene_id"]

library(ggplot2)
library(cowplot)
p.allelic.fraction <-
  ggplot(altratio.melt.ss3[!is.na(value) & detected == T & chr %in% c("16", "X")], aes(x=factor(pos, levels=unique(sort(pos))), y=reorder(variable,value,mean), fill=1-value)) +
    geom_raster() +
    facet_grid(growth+sex~chr, scales = "free", space = "free") +
    scale_fill_distiller(palette="Spectral", name="Maternal fraction") +
    theme_void()

ggsave2("plots/ss3_allelic_fraction.pdf",width = 3, height = 3, p.allelic.fraction)