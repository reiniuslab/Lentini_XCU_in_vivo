suppressMessages(source("R/load_data_smartseq3.R"))

# load escapees
x.escapees <- fread("x_escapees.tsv")

tpm.melt.ss3[expressed == T & day != 0 & x.status %in% c("XaXa","XaXi","XiXa"), fc := value / mean(value[x.status == "XaXa"],na.rm=T), by=c("l1","gene")]
tpm.melt.fcavg.ss3 <- tpm.melt.ss3[,mean(fc,na.rm=T), by=c("l1","gene","chrx","x.status")]

# p values
tpm.melt.fcavg.ss3[!gene %in% x.escapees$gene & x.status %in% c("XaXi","XiXa"), wilcox.test(V1~chrx)$p.value, by=c("l1","x.status")]

# plot
library(ggplot2)
library(cowplot)
p.fc.allele.ss3 <-
  ggplot(tpm.melt.fcavg.ss3[!gene %in% x.escapees$gene & x.status %in% c("XaXi","XiXa")], aes(x=log2(V1),fill=paste(chrx,l1=="c57"),col=paste(chrx,l1=="c57") )) +
    geom_density(alpha = 0.25) +
    facet_wrap(~x.status) +
    labs(y="Density", x="log2 fold change (vs. XaXa)") +
    coord_cartesian(xlim=c(-3,3)) +
    scale_color_brewer(palette="Paired", name="XmXp", labels=c("Autosomes:Paternal", "Autosomes:Maternal","ChrX:Paternal","ChrX:Maternal"), aesthetics = c("col","fill")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_foldchange_allele.pdf",width = 4,height = 3,p.fc.allele.ss3)