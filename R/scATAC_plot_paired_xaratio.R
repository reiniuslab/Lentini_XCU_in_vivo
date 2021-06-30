suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

## scATAC
se.gene.total <- getMatrixFromProject(ap.total.filt, "GeneScoreMatrix", logFile = NULL)
se.gene.total %<>% .[,order(colnames(se.gene.total))]

gs.all.melt <- as.data.table(melt(as.matrix(assay(se.gene.total))))
gs.all.melt[, chrom := as.character(rowData(se.gene.total)$seqnames)[Var1]]
gs.all.melt[, sex := colData(se.gene.total)[Var2,"Sex"]]
gs.all.melt[, x.status := colData(se.gene.total)[Var2,"x.status"]]

# X:A ratios
gs.xa <- gs.all.melt[, mean(value[chrom == "chrX"]) / mean(value[chrom != "chrX"]), by=c("Var2", "sex", "x.status")]

## SS3
tpm.all.melt.paired <- as.data.table(melt(tpm.all.paired))
tpm.all.melt.paired[meta.paired, c("sex", "x.status") := list(Sex, x.status), on = "Var2 == Sample_ID"]
tpm.all.melt.paired[, chrom := gene.anno.ss3[match(tpm.all.melt.paired$Var1, gene_id), paste0("chr", chrom)]]

# X:A ratios
tpm.xa.paired <- tpm.all.melt.paired[value > 1, median(value[chrom == "chrX"]) / median(value[chrom != "chrX"]), by=c("Var2", "sex", "x.status")]

## Combined
xa.comb.paired <- rbindlist(list(
  "atac" = gs.xa,
  "rna" = tpm.xa.paired
), idcol=T)
xa.comb.paired[, rel := V1 / median(V1[x.status == "XaXa"], na.rm = T), by=".id"]
xa.comb.paired[!is.na(x.status), x.status.simple := sapply(x.status, function(x) paste0(sort(c(substr(x, 1,2), substr(x, 3,4))), collapse = "") )]

library(ggplot2)
library(cowplot)
p.scatac.paired.xa <-
ggplot(xa.comb.paired[!is.na(x.status) & sex == "Female" | (sex == "Male" & x.status.simple == "XaXi")], aes(x=factor(x.status.simple, levels=c("XaXa", "XaXs", "XsXa", "XaXi", "XiXa")), y=rel, col=.id)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(0.5, 1), lty=2, col="grey") +
  facet_grid(~sex, scales="free_x", space="free_x") + 
  labs(y="Relative X:A ratio", x=NULL, col=NULL) + 
  coord_cartesian(ylim=c(0,2)) + 
  scale_color_brewer(palette="Set1", direction = -1) +
  theme_cowplot() + 
  theme(strip.background = element_blank())

ggsave2("plots/scatac_xa_ratio.pdf", width = 4, height = 4, p.scatac.paired.xa)