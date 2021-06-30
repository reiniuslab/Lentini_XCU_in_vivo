library(data.table)
# File name examples <name, allele, resolution, chromosome>: 
# D0_G1_1mb_chrX.txt
# G1 = 129, G2 = CAST
## load Hi-C contacts
fls.hic <- list.files("Froberg/dump", full.names = T)
ls.hic <- lapply(fls.hic, fread)
names(ls.hic) <- gsub("\\..*", "", basename(fls.hic))
dat.hic <- rbindlist(ls.hic, idcol = T)
colnames(dat.hic) <- c("name","x","y","value")
dat.hic[, c("name","allele","res","chrom") := tstrsplit(gsub("\\..*","",name), "_")]
dat.hic[, day := as.numeric(gsub("^D","",name))]

## correct Hi-C data for depth
dat.hic[, avgsignal := sum(value, na.rm=T)/1e6, by=c("day", "res", "chrom", "allele")]
dat.hic[, normsignal := value/avgsignal, by = c("name", "allele", "chrom", "res")]

## load pearson correlations
fls.pear <- list.files("Froberg/pearsons", full.names = T)
ls.pear <- lapply(fls.pear, function(x) as.data.table(melt(as.matrix(fread(x)))) )
names(ls.pear) <- gsub("\\..*", "", basename(fls.pear))
dat.pear <- rbindlist(ls.pear, idcol = T)
colnames(dat.pear) <- c("name","x","y","value")
dat.pear[, c("name","allele","res","chrom") := tstrsplit(gsub("\\..*","",name), "_")]
dat.pear[, day := as.numeric(gsub("^D","",name))]
dat.pear[, y := as.numeric(gsub("V","",y,fixed = T))] # fix positions

## load eigenvectors
fls.eigen <- list.files("Froberg/eigen", full.names = T)
ls.eigen <- lapply(fls.eigen, fread)
names(ls.eigen) <- gsub("\\..*", "", basename(fls.eigen))
dat.eigen <- rbindlist(ls.eigen, idcol = T)
colnames(dat.eigen) <- c("name","value")
dat.eigen[, c("name","allele","res","chrom") := tstrsplit(gsub("\\..*","",name), "_")]
dat.eigen[, day := as.numeric(gsub("^D","",name))]
dat.eigen[, x := seq_along(value)-1, by = c("name", "allele", "res", "chrom")]
dat.eigen[, chrom := factor(chrom, levels = paste0("chr", c(5,7,9,13,"X")))]

## gene density information
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
gn <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
si <- keepStandardChromosomes(SeqinfoForBSGenome(BSgenome.Mmusculus.UCSC.mm10))
bins <- tileGenome(si, tilewidth = 1e6, cut.last.tile.in.chrom = T)
score(bins) <- countOverlaps(bins, gn)

dat.gndensity <- as.data.table(bins)
dat.gndensity[,start.fix := (start - 1)/1e6]

## correlation between eigenvectors and gene density
dat.eigen[dat.gndensity, gn.density := score, on = c("chrom == seqnames", "x == start.fix")]
dat.eigen[allele == "G2" & res == "1mb", cor(value, gn.density, use = "na.or.complete"), by = c("chrom", "name")]

## rna
dat.rna <- fread("data/diffexpr_xaxi_ss3.tsv")
dat.rna[, chrom := factor(chrom, levels = c(1:19,"X"))]
dat.rna[, mbp := round(pos / 1e6)]

## PLOT
library(ggplot2)
library(cowplot)
library(scales)
## Hi-C contact maps
p.hic.chrx <-
  ggplot() +
  geom_raster(data = dat.hic[chrom == "chrX" & res == "1mb"], aes(x/1e6, y/1e6, fill=normsignal )) +
  geom_raster(data = dat.hic[chrom == "chrX" & res == "1mb"], aes(y/1e6, x/1e6, fill=normsignal )) +
  #annotate("rect", ymin=c(4,73) , xmin=c(4,73), ymax=c(73,170), xmax=c(73,170), col="black", lty=2, fill=NA) +
  facet_grid(day~allele) +
  ggtitle("Chromosome X") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "OrRd", limits=c(10,100),direction = 1, oob=squish) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.ticks.length=unit(.1, "cm"), aspect.ratio = 1 )

p.hic.chr5 <-
  ggplot() +
  geom_raster(data = dat.hic[chrom == "chr5" & res == "1mb"], aes(x/1e6, y/1e6, fill=normsignal )) +
  geom_raster(data = dat.hic[chrom == "chr5" & res == "1mb"], aes(y/1e6, x/1e6, fill=normsignal )) +
  facet_grid(day~allele) +
  ggtitle("Chromosome 5") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "OrRd", limits=c(10,100), direction = 1, oob=squish) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.ticks.length=unit(.1, "cm"), aspect.ratio = 1 )

p.hic.chr7 <-
  ggplot() +
  geom_raster(data = dat.hic[chrom == "chr7" & res == "1mb"], aes(x/1e6, y/1e6, fill=normsignal )) +
  geom_raster(data = dat.hic[chrom == "chr7" & res == "1mb"], aes(y/1e6, x/1e6, fill=normsignal )) +
  facet_grid(day~allele) +
  ggtitle("Chromosome 7") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "OrRd", limits=c(10,100), direction = 1, oob=squish) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.ticks.length=unit(.1, "cm"), aspect.ratio = 1 )

p.hic.chr9 <-
  ggplot() +
  geom_raster(data = dat.hic[chrom == "chr9" & res == "1mb"], aes(x/1e6, y/1e6, fill=normsignal )) +
  geom_raster(data = dat.hic[chrom == "chr9" & res == "1mb"], aes(y/1e6, x/1e6, fill=normsignal )) +
  facet_grid(day~allele) +
  ggtitle("Chromosome 9") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "OrRd", limits=c(10,100), direction = 1, oob=squish) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.ticks.length=unit(.1, "cm"), aspect.ratio = 1 )

p.hic.chr13 <-
  ggplot() +
  geom_raster(data = dat.hic[chrom == "chr13" & res == "1mb"], aes(x/1e6, y/1e6, fill=normsignal )) +
  geom_raster(data = dat.hic[chrom == "chr13" & res == "1mb"], aes(y/1e6, x/1e6, fill=normsignal )) +
  facet_grid(day~allele) +
  ggtitle("Chromosome 13") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "OrRd", limits=c(10,100), direction = 1, oob=squish) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.ticks.length=unit(.1, "cm"), aspect.ratio = 1 )

ggsave2("plots/hic_contact_maps.pdf", width = 22, height = 6,
plot_grid(nrow = 1, align = "hv",
  p.hic.chr5, p.hic.chr7, p.hic.chr9, p.hic.chr13, p.hic.chrx
))

## plot pearson correlations
p.pear.chrx <-
  ggplot() +
  geom_raster(data = dat.pear[chrom == "chrX" & allele == "G2" & res == "1mb" & !is.na(value)], aes(x, y, fill=value )) +
  facet_grid(day~allele, scales="free", space="free") +
  ggtitle("Chromosome X") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "RdBu", limit=c(-0.5,0.5), name = "r", na.value = NA, oob=squish) +
  theme_void() +
  theme(aspect.ratio = 1)

p.pear.chr5 <-
  ggplot() +
  geom_raster(data = dat.pear[chrom == "chr5" & allele == "G2" & res == "1mb" & !is.na(value)], aes(x, y, fill=value )) +
  facet_grid(day~allele, scales="free", space="free") +
  ggtitle("Chromosome 5") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "RdBu", limit=c(-0.5,0.5), name = "r", na.value = NA, oob=squish) +
  theme_void() +
  theme(aspect.ratio = 1)

p.pear.chr7 <-
  ggplot() +
  geom_raster(data = dat.pear[chrom == "chr7" & allele == "G2" & res == "1mb" & !is.na(value)], aes(x, y, fill=value )) +
  facet_grid(day~allele, scales="free", space="free") +
  ggtitle("Chromosome 7") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "RdBu", limit=c(-0.5,0.5), name = "r", na.value = NA, oob=squish) +
  theme_void() +
  theme(aspect.ratio = 1)

p.pear.chr9 <-
  ggplot() +
  geom_raster(data = dat.pear[chrom == "chr9" & allele == "G2" & res == "1mb" & !is.na(value)], aes(x, y, fill=value )) +
  facet_grid(day~allele, scales="free", space="free") +
  ggtitle("Chromosome 9") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "RdBu", limit=c(-0.5,0.5), name = "r", na.value = NA, oob=squish) +
  theme_void() +
  theme(aspect.ratio = 1)

p.pear.chr13 <-
  ggplot() +
  geom_raster(data = dat.pear[chrom == "chr13" & allele == "G2" & res == "1mb" & !is.na(value)], aes(x, y, fill=value )) +
  facet_grid(day~allele, scales="free", space="free") +
  ggtitle("Chromosome 13") +
  scale_y_reverse() +
  scale_fill_distiller(palette = "RdBu", limit=c(-0.5,0.5), name = "r", na.value = NA, oob=squish) +
  theme_void() +
  theme(aspect.ratio = 1)

ggsave2("plots/hic_pearson_correlation.pdf", width = 22, height = 6,
plot_grid(nrow = 1, align = "hv",
  p.pear.chr5, p.pear.chr7, p.pear.chr9, p.pear.chr13, p.pear.chrx
))

## plot eigenvectors
p.eigen.area <-
ggplot(dat.eigen[allele == "G2" & res == "1mb"], aes(y=value, x=x)) +
  geom_col() +
  facet_grid(day~chrom, scales = "free_x", space= "free_x") +
  coord_cartesian(ylim=c(-0.15,0.15)) +
  labs(y="PC1", x="Chromosome position (Mbp)") +
  theme_cowplot() +
  theme(strip.background = element_blank())

p.eigen.gndensity <-
ggplot(dat.gndensity[seqnames %in% paste0("chr", c(5,7,9,13,"X"))], aes(y=score, x=start.fix)) +
  geom_col() +
  facet_grid(0~seqnames, scales = "free", space="free") +
  coord_cartesian(ylim=c(0,60)) +
  labs(y="Genes", x="Chromosome position (Mbp)") +
  theme_cowplot() +
  theme(strip.background = element_blank())

p.eigen.rna <- 
ggplot(dat.rna[l1 == "cast" & x.status == "XiXa" & chrom %in% c(5,7,9,13,"X")], aes(x=pos/1e6, y=log2(fc) )) + 
  geom_point() + 
  geom_hline(yintercept = 0, lty=2, col="grey") +
  facet_grid(~chrom, scales = "free_x", space = "free_x") + 
  coord_cartesian(ylim=c(-3,3)) +
  labs(y="log2 fold change (CAST)", x="Chromosome position (Mbp)") +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/hic_eigen_area.pdf", height = 7, width = 10,
plot_grid(ncol=1, align="v", axis = "trbl",rel_heights = c(2,1,1),
  p.eigen.area,
  p.eigen.gndensity,
  p.eigen.rna
))