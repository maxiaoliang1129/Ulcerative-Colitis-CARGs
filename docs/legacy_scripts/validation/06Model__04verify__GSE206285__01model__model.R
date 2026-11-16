

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\06Model\\04verify\\GSE206285\\01model'
setwd(path)

#data
rt<-read.table('normalize.txt',header = T,row.names=1,sep='\t')
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')

rt<-rt[hub$V1,]
rt<-t(rt)
rt<-data.frame(rt)

rt$Sample<-rownames(rt)

rt<-merge(group,rt)
rownames(rt)<-rt[,1]
rt<-rt[,-1]
colnames(rt)[1]<-'Group'
table(rt$Group)
#control      UC 
#     18     550 

#导出验证集
exp<-rt
exp[which(exp$Group == 'control'),]$Group<-0
exp[which(exp$Group == 'UC'),]$Group<-1
write.table(exp,file="GSE206285_data.txt",sep="\t",row.names=T,quote=F)
