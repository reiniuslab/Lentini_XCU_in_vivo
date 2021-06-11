suppressMessages(source("R/load_data_smartseq3.R"))

# average data per chromosome
tpm.melt.avg.chr.ss3 <- tpm.melt.ss3[expressed == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample_id","sex","l1","day","chr","x.status")]

# average data
tpm.melt.avg.ss3 <- tpm.melt.ss3[expressed == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample_id","sex","l1","day","chrx","x.status")]

tpm.melt.all.avg.ss3 <- tpm.melt.all.ss3[expressed == T & chr != "chrY", mean(value,na.rm=T),by=c("sample_id","sex","day","chrx","x.status")]

## plot
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

# x status fraction per day
p.xstate.ss3 <-
ggplot(meta.ss3[x.status %in% c("XaXi","XiXa","XaXa")], aes(x=factor(day), fill=x.status)) +
  geom_bar(position="fill") +
  facet_grid(~sex, scales = "free_x", space = "free_x") +
  labs(x="Day", y="Fraction") +
  scale_fill_brewer(palette="Set1", name="XmXp") +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/ss3_xstate_fraction.pdf",width = 3,height = 3,p.xstate.ss3)

# tpm allele per chromosome
p.exprs.allele.chr.ss3 <- 
  ggplot(tpm.melt.avg.chr.ss3[!is.na(x.status) & x.status %in% c("XaXi","XiXa","XaXa")], aes(y=V1, x=factor(gsub("chr","",chr),level=c(1:19,"X")), col=l1=="c57" )) +
    stat_summary(fun.data="median_cl_boot") +
    facet_grid(sex+x.status~.) +
    labs(x="Chromosome", y="Expression (TPM)") +
    coord_cartesian(ylim=c(0,40)) +
    scale_color_brewer(palette="Set1", name=NULL, labels = c("Paternal", "Maternal") ) +
    theme_cowplot() +
    theme(strip.background = element_blank() )

ggsave2("plots/ss3_x_expression_chr.pdf", width = 6, height = 4, p.exprs.allele.chr.ss3)

# tpm global over days
p.exprs.all.tpm.day.ss3 <-  
  ggplot(tpm.melt.all.avg.ss3, aes(y=V1, x=factor(day), col=chrx)) +
    geom_boxplot(outlier.stroke = NA) +
    facet_grid(~sex, scales = "free_x", space = "free_x") +
    labs(x="Day", y="ChrX expression (TPM)") +
    coord_cartesian(ylim=c(0,50)) +
    scale_color_brewer(palette="Set1", name = NULL, labels = c("Autosomes", "ChrX")) +
    theme_cowplot() +
    theme(strip.background = element_blank())

ggsave2("plots/ss3_x_expression_all_tpm_day.pdf",width = 4,height = 3,p.exprs.all.tpm.day.ss3)

# tpm global over x status
p.exprs.all.tpm.xstatus.ss3 <-
  ggplot(tpm.melt.all.avg.ss3[x.status %in% c("XaXa","XaXi","XiXa") & chrx == T], aes(y=V1, x=x.status)) +
    geom_boxplot(outlier.stroke = NA, outlier.alpha = 0.1) +
    facet_grid(~sex, scales = "free_x", space = "free_x") +
    labs(x=NULL, y="Total X expression (TPM)") +
    coord_cartesian(ylim=c(0,50)) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/ss3_x_expression_all_tpm_xstatus.pdf",width = 3,height = 3,p.exprs.all.tpm.xstatus.ss3)

# xaxa expression per day
p.xexprs.all.tpm.s3 <-
  ggplot(tpm.melt.all.avg.ss3[!is.na(x.status) & x.status %in% c("XaXi","XiXa","XaXa") & chrx == T], aes(y=V1, x=factor(day), col="Total")) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    facet_grid(~sex+x.status, scales = "free_x", space = "free_x") +
    labs(x=NULL, y="Total X expression (TPM)") +
    coord_cartesian(ylim=c(0,50)) +
    theme_cowplot() +
    theme(strip.background = element_blank())

p.xexprs.tpm.s3 <-
  ggplot(tpm.melt.avg.chr.ss3[!is.na(x.status) & x.status %in% c("XaXi","XiXa","XaXa") & chr == "chrX"], aes(y=V1, x=factor(day), col=l1=="c57" )) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    facet_grid(~sex+x.status, scales="free_x", space="free_x") +
    labs(x="Day", y="Expression (TPM)") +
    coord_cartesian(ylim=c(0,50)) +
    scale_color_brewer(palette="Paired", name="ChrX", labels = c("Paternal", "Maternal") ) +
    theme_cowplot() +
    theme(strip.background = element_blank() )

ggsave2("plots/ss3_x_expression_by_day.pdf",width = 6,height = 4,
  plot_grid(align = "hv", ncol = 1,
    p.xexprs.all.tpm.s3,
    p.xexprs.tpm.s3
  )
)
