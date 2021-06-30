suppressMessages(source("R/load_data_cellrep.R"))

refratio.x.melt <- as.data.table(melt(as.matrix(counts.c57 / (counts.c57 + counts.cast))))
refratio.x.melt %<>% add.info(refseq.gene, meta)
refratio.x.melt[refseq.gene, pos := txStart, on="gene == name2"]

refratio.x.melt[,detected := mean(!is.na(value))>0.1, by="gene"]
refratio.x.melt[, maternal.fraction := ifelse(maternal == "C57", value, 1-value)]

library(ggplot2)
library(cowplot)
p.allelicfraction.x <-
  ggplot(refratio.x.melt[chr %in% c("chr16", "chrX") & !is.na(value) & sex == "F" & detected == T], aes(x=factor(pos, levels=unique(sort(pos))), y=reorder(sample,maternal.fraction,mean), fill = maternal.fraction)) +
    geom_raster() +
    scale_fill_distiller(palette="Spectral") +
    facet_grid(lineage~chr, space="free", scales="free") +
    theme_void()

ggsave2("plots/postimplantation_allelic_fraction_x.pdf", height = 3,width = 4, p.allelicfraction.x)