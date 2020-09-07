###q value FDR correction for nominal pass####
tissue=commandArgs(T)[1]
d=read.table(paste0(tissue,".nominals.txt.gz"),hea=F,stringsAsFactors=F)
colnames(d)<-c("pid","sid","dst","pval","slope")

####Storey & Tibshirani correction
library(qvalue)
library(dplyr)
dd<-
  d%>%
  group_by(pid)%>%
  mutate(pval.adj = p.adjust (pval, method='BH'))
write.table(dd[which(dd$pval.adj <= 0.05), ], paste0(tissue,".nominals.BH.FDR.txt"), quote=F, row.names=F, col.names=T)
