

rm(list=ls())

#estimate计算分值

library(limma)
library(estimate)

#path
path<-'D:\\LZ-25234\\UC\\11Immue\\02estimate'
setwd(path)

#data
rt<-read.table('normalize.txt',sep = '\t',header = T,row.names = 1)
Risk_sample<-read.table('group.txt',header = T,sep = '\t',check.names=F)

Risk_sampleid<-Risk_sample$Sample
rt<-rt[,Risk_sampleid]
dim(rt)

#处理矩阵
exp<-rt
dim(exp)

exp[1:3,1:3]

#Deduplication
counts<-exp
counts_delrep<-counts
counts_delrep<-avereps(counts_delrep)
Exp<-counts_delrep
dim(Exp)

exp = Exp
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)

write.table(data,file = "uniq.symbol.txt",sep = "\t",quote = F,col.names = T,row.names = T)

filterCommonGenes(input.f = "uniq.symbol.txt",
                  output.f = "commonGenes.gct",id = "GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",output.ds = "estimateScore.gct")

scores<-read.table("estimateScore.gct", skip = 2, header = T)

rownames(scores)<-scores[,1]
scores<-t(scores[,3:ncol(scores)])
rownames(scores)<-gsub("\\.", "\\-", rownames(scores))
out<-rbind(ID = colnames(scores), scores)

write.table(out,file = "scores.txt",sep = "\t",quote = F,col.names = F)

