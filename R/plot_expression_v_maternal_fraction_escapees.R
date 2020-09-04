suppressMessages(source("R/load_data_cellrep.R"))

# load escapee
x.escapees <- fread("x_escapees.tsv")

tpm.melt.avg.esc <- tpm.melt[expressed_lineage == T & gene %in% x.escapees$gene, base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx","maternal","c57.x.frac")]
tpm.melt.avg.esc[, maternal.x.frac := ifelse(maternal == "C57",c57.x.frac,1-c57.x.frac)]

# plot
library(ggplot2)
library(cowplot)
p.exprs.m.frac.esc <- 
  ggplot(na.omit(tpm.melt.avg.esc[sex == "F" & lineage == "EPI"]),aes(y = V1, x = maternal.x.frac, col = l1==tolower(maternal), fill = l1==tolower(maternal) )) +
  geom_vline(xintercept = 0.5, col=pal.t10[8], lty=2) +
  geom_point(stroke=NA, alpha=0.1) +
  geom_smooth(method="lm") +
  labs(x = "Maternal X-fraction", y = "Expression (TPM)") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,100)) +
  expand_limits(y=0) +
  scale_fill_brewer(palette="Paired", name = NULL, labels=c("Paternal", "Maternal") ) +
  scale_color_brewer(palette = "Paired", name = NULL, labels=c("Paternal", "Maternal")) +
  theme_cowplot() +
  theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/m_fraction_expression_escapees.pdf", width = 4, height = 4,
  p.exprs.m.frac.esc
)