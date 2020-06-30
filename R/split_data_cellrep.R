suppressMessages(source("R/load_data_cellrep.R"))

correct.values <- function(x,all){
  # corrects allelic values depending on overall detection
  # if the gene is detected => 0, else missing
  x[which(is.na(x),arr.ind = T)] <- 0
  x[which(all == 0,arr.ind = T)] <- ""
  x <- x[apply(all,1,mean) != 0,]
  return(x)
}

# correct values and export split data
tpm.allele.corr <-  lapply(tpm.allele, correct.values, all = tpm.all)

sapply(levels(spl.idx),function(x){
  write.table(tpm.allele.corr$c57[,spl.idx == x],paste0("kinetics/input/tpm/tpm.C57.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(tpm.allele.corr$cast[,spl.idx == x],paste0("kinetics/input/tpm/tpm.CAST.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})