suppressMessages(source("R/load_data_smartseq3.R"))

# average data
tpm.melt.avg.ss3 <- tpm.melt.ss3[expressed == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample_id","sex","l1","day","chrx","x.status")]

tpm.melt.all.avg.ss3 <- tpm.melt.all.ss3[expressed == T & chr != "chrY", mean(value,na.rm=T),by=c("sample_id","sex","day","chrx","x.status")]
umi.melt.all.avg.ss3 <- relumi.melt.all.ss3[expressed == T & chr != "chrY", mean(value,na.rm=T),by=c("sample_id","sex","day","chrx","x.status")]

# P-values
tpm.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & chrx == T, pairwise.wilcox.test(V1, paste(sex, x.status), "fdr")[[3]]]
umi.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & day == 0 & chrx == T, pairwise.wilcox.test(V1, paste(sex, x.status), "fdr")[[3]]]

## Plot
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
# x fraction over day
p.xfrac.ss3 <- 
ggplot(meta.ss3, aes(y=c57.x.frac, x=factor(day))) +
  geom_violin(scale = "width") +
  geom_quasirandom(alpha=0.1, stroke=NA) +
  geom_text(data=meta.ss3[!is.na(c57.x.frac), list(.N,c57.x.frac = 1), by=c("sex","day")], aes(label=N)) +
  facet_grid(~sex, scales = "free_x", space = "free_x") +
  labs(x="Day", y="Maternal X-fraction") +
  coord_cartesian(ylim=c(0,1)) +
  theme_cowplot() +
  theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_maternal_x_fraction.pdf",width = 3,height = 3,p.xfrac.ss3)

# tpm allele over x status
p.exprs.allele.tpm.ss3 <-
  ggplot(tpm.melt.avg.ss3[!is.na(x.status) & x.status %in% c("XaXi","XiXa","XaXa")], aes(y=V1, x=factor(day), col=paste(chrx, l1=="c57"))) +
    geom_boxplot(outlier.alpha = 0.1, outlier.stroke = NA) +
    facet_grid(~sex+x.status, scales = "free_x", space = "free_x") +
    labs(x="Day", y="Expression (TPM)") +
    coord_cartesian(ylim=c(0,100)) +
    scale_color_brewer(palette="Paired", name=NULL, labels=c("Autosome:Pat","Autosome:Mat","ChrX:Pat","ChrX:Mat")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression.pdf",width = 6,height = 3,p.exprs.allele.tpm.ss3)

# tpm global over x status
p.exprs.all.tpm.ss3 <-
  ggplot(tpm.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & chrx == T], aes(y=V1, x=factor(day), col=chrx)) +
    geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
    geom_hline(lty=2, col="grey",yintercept = tpm.melt.all.avg.ss3[x.status == "XaXa" & chrx == T, median(V1, na.rm = T), by=day][,median(V1)]) +
    facet_grid(~sex+x.status, scales = "free_x", space = "free_x") +
    labs(x="Day", y="Expression (TPM)") +
    coord_cartesian(ylim=c(0,100)) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_all_tpm.pdf",width = 5,height = 3,p.exprs.all.tpm.ss3)

# umi global over x status
p.exprs.all.umi.ss3 <-
  ggplot(umi.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & day == 0 & chrx == T], aes(y=V1, x=x.status, col=chrx)) +
    geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
    geom_hline(lty=2, col="grey",yintercept = umi.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & day == 0 & chrx == F, median(V1, na.rm = T)/2]) +
    facet_grid(~day+sex, scales = "free_x", space = "free_x") +
    labs(x=NULL, y="Relative UMIs") +
    coord_cartesian(ylim=c(0,1)) +
    scale_color_brewer(palette="Set1", name=NULL, labels="ChrX") +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_all_umi.pdf",width = 2,height = 3,p.exprs.all.umi.ss3)
