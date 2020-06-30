suppressMessages(source("R/load_data_cellrep.R"))

##

load.kinetics <- function(file){
  require(data.table)
  dat <- fread(file)
  names(dat)[1] <- "gene"
  dat.melt <- melt(dat, id.vars="gene")
  dat.melt[refseq.gene, chr := chrom, on = "gene == name2"] # add chromosome annotations
  dat.melt[chr != "chrY", chrx := chr == "chrX"]
  dat.melt[, c("sex", "lineage", "x.status", "allele") := tstrsplit(variable, "_", fixed=TRUE)] # split sample annotations
  dat.melt[meta[, .N , by=c("Sex", "Lineage", "x.status")], N := N, on = c("sex == Sex", "lineage == Lineage", "x.status == x.status")]
  dat.melt[, activex := F]
  dat.melt[substr(x.status,1,2) == "Xa" & allele == "C57", activex := T]
  dat.melt[substr(x.status,3,4) == "Xa" & allele == "CAST", activex := T]
  dat.melt[,value.rel := value / median(value[chrx == F], na.rm = T), by = c("sex","lineage","x.status","allele")]
  return(dat.melt)
}

##

bf.melt <- load.kinetics("kinetics/output/tpm/split/postimplantation_bf.tsv")
bs.melt <- load.kinetics("kinetics/output/tpm/split/postimplantation_bs.tsv")

library(ggplot2)
library(cowplot)
p.bf.epi <-
  ggplot(bf.melt[lineage == "EPI" & chrx == T & N > 20 & activex == T], aes(y=log2(value.rel), x=factor(paste(sex,x.status),levels=paste(rep(c("F","M"),each=5),c("XaXa","XaXs","XsXa","XaXi","XiXa"))), x2= allele )) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    stat_summary(fun.data = "median_cl_boot", position = position_dodge(0.5), stroke = NA) +
    labs(y="log2 burst frequency\n(relative to autosomes)") +
    coord_cartesian(ylim = c(-1,1)) +
    scale_color_brewer(palette = "Set1", name = "XC57XCAST") +
    facet_grid(~lineage, scales = "free_x", space = "free_x") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.bf.extra <- 
  ggplot(bf.melt[lineage != "EPI" & chrx == T & N > 20 & activex == T], aes(y=log2(value.rel), x=factor(paste(sex,x.status),levels=paste(rep(c("F","M"),each=5),c("XaXa","XaXs","XsXa","XaXi","XiXa"))), x2= allele )) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    stat_summary(fun.data = "median_cl_boot", position = position_dodge(0.5), stroke = NA) +
    labs(y="log2 burst frequency\n(relative to autosomes)") +
    coord_cartesian(ylim = c(-1,1)) +
    scale_color_brewer(palette = "Set1", name = "XC57XCAST") +
    facet_grid(~lineage, scales = "free_x", space = "free_x") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.bs.epi <-
  ggplot(bs.melt[lineage == "EPI" & chrx == T & N > 20 & activex == T], aes(y=log2(value.rel), x=factor(paste(sex,x.status),levels=paste(rep(c("F","M"),each=5),c("XaXa","XaXs","XsXa","XaXi","XiXa"))), x2= allele )) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    stat_summary(fun.data = "median_cl_boot", position = position_dodge(0.5), stroke = NA) +
    labs(y="log2 burst size\n(relative to autosomes)") +
    coord_cartesian(ylim = c(-1,1)) +
    scale_color_brewer(palette = "Set1", name = "XC57XCAST") +
    facet_grid(~lineage, scales = "free_x", space = "free_x") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.bs.extra <- 
  ggplot(bs.melt[lineage != "EPI" & chrx == T & N > 20 & activex == T], aes(y=log2(value.rel), x=factor(paste(sex,x.status),levels=paste(rep(c("F","M"),each=5),c("XaXa","XaXs","XsXa","XaXi","XiXa"))), x2= allele )) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    stat_summary(fun.data = "median_cl_boot", position = position_dodge(0.5), stroke = NA) +
    labs(y="log2 burst size\n(relative to autosomes)") +
    coord_cartesian(ylim = c(-1,1)) +
    scale_color_brewer(palette = "Set1", name = "XC57XCAST") +
    facet_grid(~lineage, scales = "free_x", space = "free_x") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

ggsave2("plots/burst_kinetics_activex.pdf", width = 4, height = 4,
  plot_grid(rel_widths = c(1.5,1), nrow = 2,
    p.bf.epi,
    p.bf.extra,
    p.bs.epi,
    p.bs.extra
  )
)