suppressMessages(source("R/load_data_smartseq3.R"))

##

load.kinetics.ss3 <- function(file){
  require(data.table)
  dat <- fread(file)
  names(dat)[1] <- "gene"
  dat.melt <- melt(dat, id.vars="gene")
  dat.melt[gene.anno.ss3, chr := paste0("chr",chrom), on = "gene == gene_id"] # add chromosome annotations
  dat.melt[chr != "chrY", chrx := chr == "chrX"]
  dat.melt[, c("sex", "growth", "x.status", "allele") := tstrsplit(variable, "_", fixed=TRUE)] # split sample annotations
  #dat.melt[meta.ss3[, .N , by=c("sex", "growth", "x.status")], N := N, on = c("sex == sex", "growth == growth", "x.status == x.status")]
  dat.melt[, activex := F]
  dat.melt[substr(x.status,1,2) == "Xa" & allele == "C57", activex := T]
  dat.melt[substr(x.status,3,4) == "Xa" & allele == "CAST", activex := T]
  dat.melt[,value.rel := value / median(value[chrx == F], na.rm = T), by = c("sex","growth","x.status","allele")]
  return(dat.melt)
}

##

tpm.bf.melt.ss3 <- load.kinetics.ss3("kinetics/output/tpm_ss3/split/tpm_ss3_bf.tsv")
tpm.bs.melt.ss3 <- load.kinetics.ss3("kinetics/output/tpm_ss3/split/tpm_ss3_bs.tsv")

umi.bf.melt.ss3 <- load.kinetics.ss3("kinetics/output/umi_ss3/split/umi_ss3_bf.tsv")
umi.bs.melt.ss3 <- load.kinetics.ss3("kinetics/output/umi_ss3/split/umi_ss3_bs.tsv")

# merge units
dat.bf.all.ss3 <- rbindlist(idcol = T,list(
  "tpm" = tpm.bf.melt.ss3,
  "umi" = umi.bf.melt.ss3
))

dat.bs.all.ss3 <- rbindlist(idcol = T,list(
  "tpm" = tpm.bs.melt.ss3,
  "umi" = umi.bs.melt.ss3
))

## p-values
dat.bf.all.ss3[activex == T & chrx ==T, pairwise.wilcox.test(value.rel, x.status, "fdr")[[3]], by=c(".id","growth")]
dat.bs.all.ss3[activex == T & chrx ==T, pairwise.wilcox.test(value.rel, x.status, "fdr")[[3]], by=c(".id","growth")]

## plot
library(ggplot2)
library(cowplot)
p.bf.ss3 <-
  ggplot(dat.bf.all.ss3[activex == T & chrx ==T], aes(y=log2(value.rel), x = x.status, col = allele, shape=.id)) +
    stat_summary(fun.data = "median_cl_boot", position=position_dodge(0.75), stroke=NA) +
    facet_grid(~sex+growth,scales = "free_x",space = "free_x") +
    labs(y="log2 burst frequency\n(relative to autosomes)", x=NULL) +
    coord_cartesian(ylim=c(-1,1)) +
    scale_color_brewer(palette="Set1", name="XmXp:", labels=c("Maternal","Paternal")) +
    scale_shape_manual(name="Unit",values=c(19,15)) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank(), axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))

p.bs.ss3 <-
  ggplot(dat.bs.all.ss3[activex == T & chrx ==T], aes(y=log2(value.rel), x = x.status, col = allele, shape=.id)) +
    stat_summary(fun.data = "median_cl_boot", position=position_dodge(0.75), stroke=NA) +
    facet_grid(~sex+growth,scales = "free_x",space = "free_x") +
    labs(y="log2 burst size\n(relative to autosomes)", x=NULL) +
    coord_cartesian(ylim=c(-1,1)) +
    scale_color_brewer(palette="Set1", name="XmXp:", labels=c("Maternal","Paternal")) +
    scale_shape_manual(name="Unit",values=c(19,15)) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank(), axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))

ggsave2("plots/ss3_kinetics_all.pdf",width = 8,height = 3,
  plot_grid(p.bf.ss3, p.bs.ss3)
)

##
p.bsbf.ss3 <-
  ggplot(cbind(dat.bf.all.ss3, value.rel.bs = dat.bs.all.ss3$value.rel)[growth=="2i" & chrx == T], aes(x=log2(value.rel), y=log2(value.rel.bs))) +
    stat_density2d(geom="polygon",aes(fill=stat(nlevel))) +
    geom_hline(yintercept = 0, col="grey", lty=2) +
    geom_vline(xintercept = 0, col="grey", lty=2) +
    facet_grid(allele~sex+x.status) +
    labs(x="log2 relative burst frequency", y="log2 relative burst size", fill="density") +
    scale_fill_viridis_c() +
    theme_cowplot() +
    theme(strip.background = element_blank())

ggsave2("plots/ss3_kinetics_relationship.pdf",width = 5,height = 4,p.bsbf.ss3)