suppressMessages(source("R/load_data_cellrep.R"))

tpm.melt.avg <- tpm.melt[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx","maternal","c57.x.frac")]
tpm.melt.avg[, maternal.x.frac := ifelse(maternal == "C57",c57.x.frac,1-c57.x.frac)]

# r2 per allele
tpm.melt.avg[sex == "F" & chrx == T & l1 == tolower(maternal) & lineage == "EPI", summary(lm(V1~maternal.x.frac))]
tpm.melt.avg[sex == "F" & chrx == T & l1 != tolower(maternal) & lineage == "EPI", summary(lm(V1~maternal.x.frac))]

# plot
library(ggplot2)
library(cowplot)
p.exprs.m.frac <- 
  ggplot(na.omit(tpm.melt.avg[sex == "F"]),aes(y = V1, x = maternal.x.frac, col = paste(chrx,l1==tolower(maternal)), fill = paste(chrx,l1==tolower(maternal)) )) +
  geom_vline(xintercept = 0.5, col=pal.t10[8], lty=2) +
  geom_point(stroke=NA, alpha=0.1) +
  geom_smooth(method="lm") +
  labs(x = "Maternal X-fraction", y = "Expression (TPM)") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,100)) +
  expand_limits(y=0) +
  facet_wrap(~lineage + chrx, ncol = 2) +
  scale_fill_brewer(palette="Paired", name = NULL, labels =c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal")) +
  scale_color_brewer(palette = "Paired", name = NULL, labels =c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal")) +
  theme_cowplot() +
  theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/m_fraction_expression.pdf", width = 4, height = 6,
  p.exprs.m.frac
)

# gene-wise correlations
tpm.melt[, maternal.x.frac := ifelse(maternal == "C57",c57.x.frac,1-c57.x.frac)]
tpm.melt[, mat := l1 == tolower(maternal)]

x.escapees <- fread("x_escapees.tsv")

tpm.melt.cor.matx <- tpm.melt[expressed_lineage == T & sex == "F" & lineage == "EPI" & !is.na(value), cor(value, maternal.x.frac),by = c("gene","mat","chrx")]
tpm.melt.cor.matx[, x.esc := gene %in% x.escapees$gene]

# ks test for escapees vs non-escapees
tpm.melt.cor.matx[chrx == T, ks.test(V1[x.esc == T],V1[x.esc == F])$p.value, by=mat]

# plot
p.cor.m.frac <-
  ggplot(tpm.melt.cor.matx[!is.na(chrx)], aes(x = V1, col= paste(chrx,mat), fill= paste(chrx,mat) )) + 
    geom_density(alpha=0.25) +
    geom_rug(data = tpm.melt.cor.matx[x.esc==T], sides = "t") +
    labs(x = "Correlation (r)") +
    scale_color_brewer(palette="Paired", name=NULL, labels=c("Autosome:Paternal","Autosome:Maternal","ChrX:Paternal","ChrX:Maternal"), aesthetics = c("col","fill")) +
    theme_cowplot() +
    guides(fill = guide_legend(nrow=2,byrow=TRUE), col = guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.position = "top")

ggsave2("plots/m_fraction_correlation.pdf",width = 4,height = 4,
  p.cor.m.frac
)