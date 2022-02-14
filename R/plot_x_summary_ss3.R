suppressMessages(source("R/load_data_smartseq3.R"))

# average data per chromosome
tpm.melt.avg.chr.ss3 <- tpm.melt.ss3[expressed == T & chr != "chrY", base::mean(value,na.rm=T, trim = 0.2),by=c("sample_id","sex","l1","day","chr","x.status")]

library(ggplot2)
library(cowplot)
p.xtotal.dodge <- 
ggplot(tpm.melt.avg.chr.ss3[chr == "chrX" & sex == "F" & x.status %in% c("XaXa", "XaXi", "XiXa")], aes(x=x.status, y=V1, fill=l1)) + 
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9), aes(col=l1)) +
  stat_summary(fun.y="mean", geom="col", position=position_dodge(0.9)) + 
  stat_summary(fun.data="mean_cl_boot", geom="linerange",position=position_dodge(0.9), aes(group=l1)) + 
  labs(y="Expression (TPM)", x=NULL) +
  coord_cartesian(ylim=c(0,45)) +
  scale_fill_brewer(palette="Paired", direction = -1, name="ChrX", aesthetics = c("col", "fill")) +
  theme_cowplot()

p.xtotal.stack <- 
ggplot(tpm.melt.avg.chr.ss3[chr == "chrX" & sex == "F" & x.status %in% c("XaXa", "XaXi", "XiXa")], aes(x=x.status, y=V1, fill=l1)) + 
  stat_summary(fun.y="mean", geom="col", position="stack", width=0.5) + 
  labs(y="Expression (TPM)", x=NULL) +
  coord_cartesian(ylim=c(0,45)) +
  scale_fill_brewer(palette="Paired", direction = -1, name="ChrX") +
  theme_cowplot()

ggsave2("plots/ss3_x_summary.pdf", height = 4, width = 4,
plot_grid(ncol = 1,
  p.xtotal.dodge,
  p.xtotal.stack
))
