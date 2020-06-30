suppressMessages(source("R/load_data_smartseq3.R"))

##

correct.values <- function(x,all){
  # corrects allelic values depending on overall detection
  # if the gene is detected => 0, else missing
  x[which(is.na(x),arr.ind = T)] <- 0
  x[which(all == 0,arr.ind = T)] <- ""
  x <- x[apply(all,1,mean) != 0,]
  return(x)
}

##

# recast to have matching rows & columns
## tpm
tpm.all.ss3 <- data.frame(dcast(tpm.melt.all.ss3, var1 ~ sample_id),row.names = 1)

tpm.allele.ss3 <- list(
  "c57" = data.frame(dcast(tpm.melt.ss3[l1 == "c57"],var1 ~ sample_id, value.var = "value"),row.names = 1),
  "cast" = data.frame(dcast(tpm.melt.ss3[l1 == "cast"],var1 ~ sample_id, value.var = "value"),row.names = 1)
)
## umi
umi.all.ss3 <- data.frame(dcast(relumi.melt.all.ss3, var1 ~ sample_id),row.names = 1)

umi.allele.ss3 <- list(
  "c57" = data.frame(dcast(umi.melt.ss3[l1 == "c57"],var1 ~ sample_id, value.var = "umi_rel"),row.names = 1),
  "cast" = data.frame(dcast(umi.melt.ss3[l1 == "cast"],var1 ~ sample_id, value.var = "umi_rel"),row.names = 1)
)

# correct values and export split data
tpm.allele.corr.ss3 <-  lapply(tpm.allele.ss3, correct.values, all = tpm.all.ss3[row.names(tpm.allele.ss3$c57), colnames(tpm.allele.ss3$c57)])
umi.allele.corr.ss3 <-  lapply(umi.allele.ss3, correct.values, all = umi.all.ss3[row.names(umi.allele.ss3$c57), colnames(umi.allele.ss3$c57)])

spl.idx.tpm.ss3 <- factor(meta.ss3[match(colnames(tpm.allele.corr.ss3$c57),sample_id), paste(sex,ifelse(day==0,"2i","Diff"),x.status,sep=".")])
spl.idx.umi.ss3 <- factor(meta.ss3[match(colnames(umi.allele.corr.ss3$c57),sample_id), paste(sex,ifelse(day==0,"2i","Diff"),x.status,sep=".")])

sapply(levels(spl.idx.tpm.ss3)[table(spl.idx.tpm.ss3)>20],function(x){
  write.table(tpm.allele.corr.ss3$c57[,spl.idx.tpm.ss3 == x],paste0("kinetics/input/tpm_ss3/tpm.C57.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(tpm.allele.corr.ss3$cast[,spl.idx.tpm.ss3 == x],paste0("kinetics/input/tpm_ss3/tpm.CAST.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})

sapply(levels(spl.idx.umi.ss3)[table(spl.idx.umi.ss3)>20],function(x){
  write.table(umi.allele.corr.ss3$c57[,spl.idx.umi.ss3 == x],paste0("kinetics/input/umi_ss3/umi.C57.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(umi.allele.corr.ss3$cast[,spl.idx.umi.ss3 == x],paste0("kinetics/input/umi_ss3/umi.CAST.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})