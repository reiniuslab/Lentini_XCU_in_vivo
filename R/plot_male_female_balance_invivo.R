suppressMessages(source("R/load_data_cellrep.R"))

# add pseudotime and active X information
pseudo <- fread("data/pseudotime.tsv")

tpm.melt[pseudo, pseudotime := Pseudotime, on="sample==LibraryName"]
tpm.melt[, l1_x := ifelse(l1 == tolower(maternal),substr(x.status.maternal,1,2),substr(x.status.maternal,3,4))]
tpm.melt[, activex := ifelse(l1_x == "Xa",T,F)]

# plot
library(ggplot2)
library(cowplot)
# over pseudotime
p.pseudo.xstatus <-
  ggplot(unique(tpm.melt[lineage=="EPI" & sex == "F",list(sample,x.status.maternal,pseudotime)]), aes(x=pseudotime,fill=x.status.maternal,)) +
    geom_area(position="fill", stat="bin", bins=21) +
    scale_fill_brewer(palette="Set1", name="XMXP") +
    theme_cowplot() +
    theme(legend.position = "top")

p.pseudo.exprs <-
  ggplot(tpm.melt[expressed_lineage == T & lineage == "EPI" & activex == T & !is.na(chrx), base::mean(value,na.rm=T,trim=0.2), by=c("sample","sex","pseudotime","l1","chrx")], aes(y = V1, x=pseudotime, col = sex, fill = sex)) +
    geom_point(stroke=NA) +
    geom_smooth(method="lm",se = F) +
    geom_smooth(method="lm",fullrange=T, lty=3) +
    coord_cartesian(ylim=c(0,80)) +
    expand_limits(y=0) +
    facet_wrap(~chrx, ncol = 1) +
    labs(y="Expression (TPM)") +
    scale_colour_brewer(palette="Set1",aesthetics = c("colour","fill")) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/mf_timeseries_pseudotime.pdf",width = 4, height = 6,
  plot_grid(rel_heights = c(0.5,1),ncol = 1,
    p.pseudo.xstatus,
    p.pseudo.exprs
  )
  
)

# over embryonic day
p.eday.xstatus <-
  ggplot(unique(tpm.melt[lineage=="EPI" & sex == "F",list(sample,x.status.maternal,embryonicday)]), aes(x=embryonicday,fill=x.status.maternal)) +
    geom_bar(position="fill") +
    scale_fill_brewer(palette="Set1", name="XMXP") +
    scale_x_continuous(breaks=c(5.5,6.0,6.5)) +
    theme_cowplot() +
    theme(legend.position = "top")

p.eday.exprs <-
  ggplot(tpm.melt[expressed_lineage == T & lineage == "EPI" & activex == T & !is.na(chrx), base::mean(value,na.rm=T,trim=0.2), by=c("sample","sex","embryonicday","l1","chrx")], aes(y = V1, x=embryonicday, col = sex, fill = sex)) +
    geom_jitter(width=0.1, stroke=NA) +
    geom_smooth(method="lm",se = F) +
    geom_smooth(method="lm",fullrange=T, lty=3) +
    coord_cartesian(ylim=c(0,80)) +
    expand_limits(y=0) +
    facet_wrap(~chrx, ncol = 1) +
    labs(y="Expression (TPM)") +
    scale_colour_brewer(palette="Set1",aesthetics = c("colour","fill")) +
    scale_x_continuous(breaks=c(5.5,6.0,6.5)) +
    theme_cowplot() +
    theme(legend.position = "top", strip.background = element_blank())

ggsave2("plots/mf_timeseries_eday.pdf",width = 4, height = 6,
  plot_grid(rel_heights = c(0.5,1),ncol = 1,
    p.eday.xstatus,
    p.eday.exprs
  )      
)