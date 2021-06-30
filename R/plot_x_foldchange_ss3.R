suppressMessages(source("R/load_data_smartseq3.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

## calculate fold changes
tpm.melt.ss3[expressed == T & day != 0 & x.status %in% c("XaXa","XaXi","XiXa"), fc := value / mean(value[x.status == "XaXa"],na.rm=T), by=c("l1","gene")]
tpm.melt.fcavg.ss3 <- tpm.melt.ss3[,list(fc = mean(fc,na.rm=T)), by=c("l1","gene","chrx","x.status")]

# tpm.melt.fcavg.ss3[gene.anno.ss3, c("pos", "chrom") := list(start, chrom), on="gene == gene_name"]
# fwrite(tpm.melt.fcavg.ss3, "data/diffexpr_xaxi_ss3.tsv", quote = F, sep = "\t")

# p values
tpm.melt.fcavg.ss3[!gene %in% x.escapees$gene & x.status %in% c("XaXi","XiXa"), wilcox.test(fc~chrx)$p.value, by=c("l1","x.status")]

## calculate average expression
tpm.melt.gnavg.ss3 <- tpm.melt.ss3[expressed == T & day != 0 & x.status %in% c("XaXa","XaXi","XiXa"), mean(value, na.rm=T), by=c("gene", "chrx", "x.status", "l1")]
tpm.cast.gnavg.ss3 <- dcast(tpm.melt.gnavg.ss3, gene+chrx+l1~x.status, value.var = "V1")

## plot
library(ggplot2)
library(cowplot)
# scatter
p.avg.c57 <-
  ggplot(tpm.cast.gnavg.ss3[l1 == "c57"], aes(x=log2(XaXa), y=log2(XaXi), col=chrx)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0 ) +
    ggtitle("C57") +
    coord_cartesian(xlim=c(-2,12), ylim=c(-2,12)) +
    scale_color_brewer(palette="Set1", name=NULL, labels=c("Autosomes","ChrX") ) +
    theme_cowplot() +
    theme(aspect.ratio = 1)

p.avg.cast <-
  ggplot(tpm.cast.gnavg.ss3[l1 == "cast"], aes(x=log2(XaXa), y=log2(XiXa), col=chrx)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0 ) +
    ggtitle("CAST") +  
    coord_cartesian(xlim=c(-2,12), ylim=c(-2,12)) +
    scale_color_brewer(palette="Set1", name=NULL, labels=c("Autosomes","ChrX") ) +
    theme_cowplot() +
    theme(aspect.ratio = 1)

ggsave2("plots/ss3_scatter_allele_expression.pdf", width = 8, height = 3,
  plot_grid(nrow = 1,
    p.avg.c57, p.avg.cast
  )
)

# density
p.dens.fc <- 
ggplot(tpm.melt.fcavg.ss3[x.status %in% c("XaXi", "XiXa")], aes(x=log2(fc), col=chrx, fill=chrx )) + 
  geom_density() + 
  geom_vline(xintercept = 0) +
  facet_grid(l1 ~ x.status) +
  coord_cartesian(xlim=c(-4, 4)) +
  scale_color_brewer(palette="Set1", aesthetics = c("colour", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank() )

ggsave2("plots/ss3_scatter_allele_fc.pdf", width = 8, height = 3, p.dens.fc)