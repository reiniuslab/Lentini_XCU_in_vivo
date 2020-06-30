#!/usr/bin/Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(TRUE)

in_dir <- args[1]
out_prefix <- args[2]

fls <- list.files(in_dir,"*_ML.tsv", full.names = T)

dat <- lapply(fls,fread)
names(dat) <- lapply(strsplit(gsub("tpm.|_.*","",basename(fls)),"\\."), function(x) paste0(x[c(2:4,1)], collapse="_"))

dat.melt <- rbindlist(dat,idcol = "sample")
names(dat.melt)[2] <- "gene"

dat.melt[,burst_size := k_syn / k_off]

bf.cast <- dcast(dat.melt[`1` == T], gene ~ sample, value.var = "k_on", fill = NA)
bs.cast <- dcast(dat.melt[`1` == T], gene ~ sample, value.var = "burst_size", fill = NA)

fwrite(bf.cast, paste0(out_prefix,"_bf.tsv"), quote = F,sep = "\t")
fwrite(bs.cast, paste0(out_prefix,"_bs.tsv"), quote = F,sep = "\t")
