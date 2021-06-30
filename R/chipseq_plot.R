library(data.table)

##

average.deeptools <- function(x){
  mat <- as.matrix(x[,-1:-6])
  mat[which(is.na(mat), arr.ind = T)] <- 0
  mat.avg <- (mat[,1:1000]+mat[,1001:2000]) / 2
  row.names(mat.avg) <- gsub(" ", "", apply(x[,1:3],1,paste0,collapse=":"))
  mat.avg.filt <- mat.avg[rowMeans(mat.avg, na.rm = T) > 0,]
  return(mat.avg.filt)
}

multiOverlap <- function(x){
  require(GenomicRanges)
  sapply(seq_along(x), function(i) mean(x[[i]] %over% Reduce(c, x[-i]) ) )
}

##

# Peak numbers
library(rtracklayer)
#fls.chip <- list.files("peaks/consensus", ".broadPeak$", full.names = T)
fls.chip <- list.files("ChIPseq/peaks/filtered", ".broadPeak$", full.names = T)

ls.chip <- lapply(fls.chip, import)
names(ls.chip) <- gsub("_peaks_.*", "", basename(fls.chip))

dat.chip <- rbindlist(lapply(ls.chip, as.data.table), idcol = T)
dat.chip[, c("data", "mod", "celline", "genotype", "timepoint", "treatment", "replicate", "allele") := tstrsplit(.id, "_")]
dat.chip[, time := as.numeric(gsub("h","",timepoint))]
dat.chip[, chrx := seqnames == "X"]

dat.chip.cnt <- dat.chip[, .N, by=c("seqnames", "allele", "time", "replicate", "mod")]
dat.chip.cnt[, N_rel := N / N[time == 0], by = c("seqnames", "allele", "replicate", "mod")]

# Overlap
dat.anno <- as.data.table(do.call(rbind, strsplit(names(ls.chip), "_")))
colnames(dat.anno) <- c("assay", "mod", "celline", "genotype", "time", "treatment", "rep", "allele")

ls.chip.red <- tapply(ls.chip, dat.anno[, paste(mod, time, allele, sep="_")], function(x) keepSeqlevels(reduce(Reduce(c, x), min.gapwidth = 500), "X", pruning.mode = "coarse") )

dat.ol <- as.data.table(melt(do.call(cbind, tapply(ls.chip.red, gsub("_.*_", "_", names(ls.chip.red)), multiOverlap))))
dat.ol[, c("mod", "allele") := tstrsplit(Var2, "_")]
dat.ol[, time := c("0h", "12h", "24h", "4h", "8h")[Var1] ]

# TSS enrichment
fls.chip.tssenirch <- list.files("ChIPseq/deeptools", "TSS_enrichment.mat.gz$", full.names = T)
ls.tssenrich.avg <- lapply(fls.chip.tssenirch, function(x){average.deeptools(fread(x, skip=1))})
names(ls.tssenrich.avg) <- gsub("\\..*","",basename(fls.chip.tssenirch))

dat.tssenrich <- as.data.table(melt(lapply(ls.tssenrich.avg, colMeans)))
dat.tssenrich[, pos := seq(-5e3, 5e3, length.out = .N), by="L1"]
dat.tssenrich[, c("mod", "time", "allele", "target", "type") := tstrsplit(L1, "_")]
dat.tssenrich[, rel := value / value[mod == "Input"], by=c("pos", "allele", "time")]

# Enhancer enrichment
fls.chip.enhenirch <- list.files("ChIPseq/deeptools", "Enhancer_enrichment.mat.gz$", full.names = T)
ls.enhenrich.avg <- lapply(fls.chip.enhenirch, function(x){average.deeptools(fread(x, skip=1))})
names(ls.enhenrich.avg) <- gsub("\\..*","",basename(fls.chip.enhenirch))

dat.enhenrich <- as.data.table(melt(lapply(ls.enhenrich.avg, colMeans)))
dat.enhenrich[, pos := seq(-5e3, 5e3, length.out = .N), by="L1"]
dat.enhenrich[, c("mod", "time", "allele", "target", "type") := tstrsplit(L1, "_")]
dat.enhenrich[, rel := value / value[mod == "Input"], by=c("pos", "allele", "time")]

## PLOT
# peaks
library(ggplot2)
library(cowplot)
p.chip.relpeaks <-
ggplot(dat.chip.cnt, aes(x=time, y=log2(N_rel), col=paste(seqnames == "X", allele == "C57"), group=paste(mod, seqnames == "X", replicate, allele) )) +
  stat_summary(fun.y="median", geom="line") +
  stat_summary(fun.y="median", geom="point", stroke=NA) +
  facet_grid(mod~replicate) +
  labs(x="Time", y="log2 relative peaks") +
  scale_x_continuous(breaks= seq(0,24,12),labels = seq(0,24,12) ) +
  coord_cartesian(ylim=c(-3, 3)) +
  scale_color_brewer(palette="Paired", name=NULL, labels=c("Autosomes:CAST", "Autosomes:C57", "ChrX:CAST", "ChrX:C57")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = "black"))

ggsave2("plots/chipseq_relativepeaks.pdf", width = 4, height = 8, p.chip.relpeaks)

# overlap
library(ggbeeswarm)
p.chip.overlap <-
ggplot(dat.ol, aes(x=factor(time, levels=c("0h","4h","8h","12h","24h")), y=value, fill=allele)) + 
  geom_col(position = "dodge") +
  facet_grid(mod~1) +
  labs(x="Time", y="Peak overlap") +
  scale_fill_brewer(palette="Set1", name="Allele") +
  expand_limits(y=0) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = "black"))

ggsave2("plots/chipseq_peakoverlap.pdf", width = 3, height = 8, p.chip.overlap)
  
# tss lineplot
p.tss.line <- 
ggplot(dat.tssenrich[mod != "Input"], aes(x=pos/1e3, y=rel, col=factor(time, levels=c("0h","4h","8h","12h","24h")), group=L1 )) +
  geom_line() +
  facet_grid(mod~allele, scales="free_y") +
  labs(y="Relative enirchment", x="Distance to TSS (kb)", col="Time (h)") +
  scale_color_viridis_d() +
  expand_limits(y=c(0,5)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/chipseq_lineplot_tss.pdf", width = 4, height = 8, p.tss.line)

# enhancer lineplot
p.enh.line <- 
ggplot(dat.enhenrich[mod != "Input"], aes(x=pos/1e3, y=rel, col=factor(time, levels=c("0h","4h","8h","12h","24h")), group=L1 )) +
  geom_line() +
  facet_grid(mod~allele, scales="free_y") +
  labs(y="Relative enirchment", x="Distance to enhancer (kb)", col="Time (h)") +
  scale_color_viridis_d() +
  expand_limits(y=c(0,5)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/chipseq_lineplot_enhancer.pdf", width = 4, height = 8, p.enh.line)
