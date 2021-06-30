suppressMessages(source("R/load_data_scATAC_scRNA_paired.R"))

##

groupXalleles <- function(ref, alt, prob = 1, x.frac ){
  if(missing(x.frac)){
    message("Calculating allelic ratio")
    total <- ref+alt
    x.avg <- rowMeans(total, na.rm=T)
    idx <- x.avg < quantile(x.avg[x.avg > 0], probs = prob, na.rm=T)
    x.frac <- colSums(ref[idx,], na.rm = T)/colSums(total[idx,], na.rm = T)
  }
  
  xci.comp <- abs(0.5-x.frac)/0.5
  
  xi <- xa <- matrix(nrow = nrow(ref), ncol = ncol(ref))
  row.names(xi) <- row.names(xa) <- row.names(ref)
  
  for(i in 1:ncol(ref)){
    if(is.na(x.frac[i])){
      xi[,i] <- NA
      xa[,i] <- NA
    }else{
      if(x.frac[i] > 0.5){
        xi[,i] <- alt[,i]
        xa[,i] <- ref[,i]
      }else{
        xi[,i] <- ref[,i]
        xa[,i] <- alt[,i]
      }
    }
  }
  
  return(list(xa = xa, xi = xi, xci = xci.comp))
}  

match.intersect <- function(x, y){
  xy <- sort(intersect(x, y))
  match(xy, x)
}

##

## Gene positions for mm10
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
gn <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
dat.gn <- as.data.table(gn)[seqnames == "chrX", .N, by=floor(start/1e6)]
dat.gn2 <- rbind(dat.gn[, list(floor, N, L1 = "xa")], dat.gn[, list(floor, N, L1 = "xi")])

## atac
# use 500 bp tiles for chrX accessibility
se.tile <- getMatrixFromProject(ap.allele.filt, "TileMatrix", binarize = T, useSeqnames = "chrX", logFile=NULL)
se.tile %<>% .[,order(colnames(se.tile))]

# group data by X alleles, use rowMeans as TileMatrix data is binary
x.atac <- groupXalleles(
  assay(se.tile)[,as.logical(se.tile$Sample == "scATAC_C57" & se.tile$Sex == "Female") & !is.na(se.tile$c57.x.frac)],
  assay(se.tile)[,as.logical(se.tile$Sample == "scATAC_CAST" & se.tile$Sex == "Female") & !is.na(se.tile$c57.x.frac)]
)
x.atac.avg <- lapply(x.atac[1:2], getGroupedMatrix, cut(x.atac$xci, c(0,0.2,0.8,1) ), fun = rowMeans, na.rm=T )

x.atac.melt <- as.data.table(melt(x.atac.avg))
x.atac.melt[, pos := rowData(se.tile)$start[Var1]] # add chromosome position
x.atac.melt[, value.norm := (value / sum(value, na.rm=T)) * 1e4] # scale to 1e4

# rolling average
x.atac.rollavg <- x.atac.melt[, list(xi = frollmean(value.norm, 5e3, na.rm=T), pos = frollmean(pos/1e6, 5e3, na.rm=T)), by= c("Var2", "L1")]
x.atac.rollavg[, mbp := round(pos)]

## rna
x.rna <- groupXalleles(
  tpm.allele.paired$c57[gene.anno.ss3[chrom == "X", intersect(row.names(tpm.allele.paired$c57), unique(gene_id))], meta.paired[Sex == "Female" & !is.na(c57.x.frac) & exclude == F, Sample_ID] ],
  tpm.allele.paired$cast[gene.anno.ss3[chrom == "X", intersect(row.names(tpm.allele.paired$cast), unique(gene_id))], meta.paired[Sex == "Female" & !is.na(c57.x.frac) & exclude == F, Sample_ID] ]
)
x.rna.avg <- lapply(x.rna[1:2], getGroupedMatrix, cut(x.rna$xci, c(0,0.2,0.8,1) ), fun = rowMeans, na.rm=T )

x.rna.melt <- as.data.table(melt(x.rna))
x.rna.melt[, pos := gene.anno.ss3[match(Var1, gene_id), start] ]
x.rna.melt[, outlier := value > median(value, na.rm = T) + mad(value,na.rm = T)*3, by=c("Var2", "L1")]
x.rna.melt %<>% .[order(pos)]

# divide into blocks to avoid averaging over large areas lacking genes
x.rna.melt[outlier == F, block := findInterval(seq_along(pos),which(c(F,diff(pos) > 9.25e6))), by=c("Var2", "L1")]
x.rna.rollavg <- x.rna.melt[outlier == F, list(xi = frollmean(value, 20, na.rm=T), pos = frollmean(pos/1e6, 20, na.rm=T)), by= c("Var2", "L1", "block")]
x.rna.rollavg[, mbp := round(pos)]

## Zylics et al. 2019 native ChIP-seq data
library(rtracklayer)
fls.chip <- list.files("Zylicz/bigwig", pattern = "ChIPseq_H3K27me3_.*_Rep1_.*.bw|ChIPseq_H3K4me1_.*_Rep1_.*.bw", full.names = T)
ls.chip <- lapply(fls.chip, import, selection = BigWigSelection(ranges = GRanges("X", IRanges(1,171031299)))) # load only chrX
names(ls.chip) <- gsub("ChIPseq_|_TX1072_WT|_Dox_Rep1|\\.bw","", basename(fls.chip))

dat.chip <- rbindlist(lapply(ls.chip, as.data.table), idcol = T)
dat.chip[, c("mod", "time", "allele") := tstrsplit(.id, "_")]

dat.chip.rollavg <- dat.chip[, list(score=frollmean(score, 2e4, na.rm=T), pos=frollmean(start/1e6, 2e4, na.rm=T)), by=c("mod","time", "allele")]
dat.chip.rollavg[, time := factor(time, levels=c("0h", "4h", "8h", "12h", "24h"))]
dat.chip.rollavg[, L1 := ifelse(allele == "C57", "xi", "xa")]
dat.chip.rollavg[, mbp := round(pos)]

## Paired XCI completion
xci.comp.paired <- 
  data.table(
    "atac" = x.atac$xci[match.intersect(gsub(".*#|_S[[:digit:]]{1,3}|D_", "", names(x.atac$xci)), gsub("R_", "", names(x.rna$xci)) )],
    "rna" = x.rna$xci[match.intersect(gsub("R_", "", names(x.rna$xci)), gsub(".*#|_S[[:digit:]]{1,3}|D_", "", names(x.atac$xci)) )],
    "rn" = names(x.rna$xci[match.intersect(gsub("R_", "", names(x.rna$xci)), gsub(".*#|_S[[:digit:]]{1,3}|D_", "", names(x.atac$xci)) )])
  )

## PLOT
library(ggplot2)
library(cowplot)
# collapse to mbp resolution for exporting
p.silencing.atac <-
ggplot(x.atac.rollavg[, mean(xi, na.rm=T), by=c("Var2", "L1", "mbp")], aes(x=mbp, y=V1, col=Var2 )) + 
  geom_line() +
  facet_grid(~L1) +
  coord_cartesian(xlim=c(5,170)) +
  scale_color_brewer(palette = "Set1", name = "XCI completion", labels = c("0-20%", "20-80%", "80-100%")) +
  labs(y="Accessibility") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"), axis.title.x = element_blank(), axis.text.x = element_blank())

p.silencing.rna <- 
ggplot(x.rna.rollavg[, mean(xi, na.rm=T), by=c("Var2", "L1", "mbp", "block")], aes(x=mbp, y=V1, col=factor(Var2) )) + 
  geom_line(lty=2, aes(group=Var2)) +
  geom_line(aes(group= paste(Var2, block) )) +
  facet_grid(~L1) +
  coord_cartesian(xlim=c(5,170)) +
  scale_color_brewer(palette = "Set1", name = "XCI completion", labels = c("0-20%", "20-80%", "80-100%")) +
  labs(y="Expression (TPM)", x="ChrX (Mbp)") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"), axis.title.x = element_blank(), axis.text.x = element_blank())

p.silencing.h3k27me3 <-
ggplot(dat.chip.rollavg[mod == "H3K27me3", mean(score, na.rm=T), by=c("time", "mbp", "L1")], aes(x=mbp, y=V1, col=time)) + 
  geom_line(aes(group=time)) +
  facet_grid(~L1) +
  labs(y="H3K27me3") +
  coord_cartesian(xlim=c(5,170), ylim=c(1,10)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"), axis.title.x = element_blank(), axis.text.x = element_blank())

p.silencing.h3k4me1 <-
ggplot(dat.chip.rollavg[mod == "H3K4me1", mean(score, na.rm=T), by=c("time", "mbp", "L1")], aes(x=mbp, y=V1, col=time)) + 
  geom_line(aes(group=time)) +
  facet_grid(~L1) +
  labs(y="H3K4me1") +
  coord_cartesian(xlim=c(5,170), ylim=c(1,10)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"), axis.title.x = element_blank(), axis.text.x = element_blank())

p.silencing.gn <- 
ggplot(dat.gn2, aes(x=floor, y=N)) +
  geom_col() +
  facet_grid(~L1) +
  annotate("rect", xmin = 103.458373, xmax=103.485233, ymin = 20, ymax=30, col="red") + # Xist location
  labs(y="Genes", x="ChrX (Mbp)") +
  coord_cartesian(xlim=c(5,170)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"))

ggsave2("plots/scatac_spatial_chrx.pdf", width = 8, height = 8,
  plot_grid(p.silencing.atac, p.silencing.rna, p.silencing.h3k4me1, p.silencing.h3k27me3, p.silencing.gn, ncol=1, align = "v", axis = "trbl", rel_heights = c(1,1,1,1,0.75))
)

##

p.xci.comp.paired <-
ggplot(melt(xci.comp.paired[rn %in% meta.paired[Sex == "Female", Sample_ID]], id.vars = "rn"), aes(y=reorder(rn,-value), x=variable, fill=value)) + 
  geom_tile() + 
  labs(y="Cells", fill="XCI completion") +
  scale_fill_viridis_c(option = "E") + 
  theme_void() + 
  theme(axis.text.x = element_text(), axis.title.y = element_text(angle = 90))

ggsave2("plots/scatac_xci_completion_paired.pdf", width=2, height = 8, p.xci.comp.paired)