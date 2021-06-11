suppressMessages(source("R/load_data_science_ss1.R"))
suppressMessages(source("R/load_data_science.R"))

##

median_cl_boot <- function(x, conf = 0.95) {
  # function from http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}

##

# average
tpm.melt.early.avg <- tpm.melt.early[expressed_lineage == T & chr != "chrY", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","lineage","sex","l1","chrx","maternal")]

tpm.melt.blast.avg <- tpm.melt.blast[expressed_lineage == T & chr != "chrY", base::mean(value,na.rm=T, trim=0.2),by=c("sample","lineage","sex","l1","chrx", "maternal")]
tpm.melt.blast.avg[,lineage := factor(lineage,levels=c("earlyblast","midblast","lateblast"))]

## plot
library(ggplot2)
library(cowplot)
p.x.expr.early <- 
  ggplot(tpm.melt.early.avg, aes(y = V1, x= lineage, x2= l1, col=paste(chrx, tolower(maternal) == l1) )) +
    stat_summary(fun.y="median", geom="line",aes(group=paste(chrx, tolower(maternal) == l1) )) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    coord_cartesian(ylim=c(0,75)) +
    ylab("Expression (TPM)") +
    facet_wrap(~sex) +
    scale_color_brewer(name = "C57xCAST:",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())

p.x.expr.blast <- 
  ggplot(tpm.melt.blast.avg,aes(y=V1,x=lineage,x2=factor(l1,levels=c("cast","c57")),col=paste(chrx, tolower(maternal) == l1))) +
    stat_summary(fun.y="median", geom="line",aes(group=paste(chrx, tolower(maternal) == l1) )) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    expand_limits(y=0) +
    ylab("Expression (TPM)") +
    coord_cartesian(ylim=c(0,75)) +
    facet_wrap(~sex) +
    scale_color_brewer(name = "C57xCAST",labels = c("Autosome:Paternal","Autosome:Maternal","chrX:Paternal","chrX:Maternal"), palette = "Paired") +
    theme_cowplot() +
    theme(legend.position = "top", axis.title.x = element_blank(), strip.background = element_blank())


ggsave2("plots/x_expression_preimplantation.pdf",height = 4, width = 18,
  plot_grid(nrow = 1, rel_widths = c(3,1),
            p.x.expr.early, p.x.expr.blast
  )
)