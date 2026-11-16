

rm(list=ls())

#path
path<-'E:\\LZ-25007\\colitis\\06verify\\GSE128682\\00data'
setwd(path)

library(limma)
library(edgeR)
library(ggvenn)
library(patchwork)
library(data.table)

#data
data<-read.table("GSE128682_DESeq2_normalized_counts.tab",header = T,row.names=1,sep = "\t")
dim(data)
#[1] 26475    44

#只要UC样本
mRNA_rt<-data[,15:ncol(data)]
mRNA_rt<-mRNA_rt[which(rowSums(mRNA_rt) > 0),]
dim(mRNA_rt)
#[1] 23777    30

#处理矩阵
merge_data<-mRNA_rt
exp<-merge_data
dim(exp)  #[1] 23777    30

exp[1:3,1:3]

#数据转换
exp<-log2(exp+1)

#Deduplication
counts<-exp
counts_delrep<-counts
counts_delrep<-avereps(counts_delrep)
Exp<-counts_delrep
dim(Exp)
#[1] 23777    30
Exp<-na.omit(Exp)
dim(Exp)

write.table(Exp,"deduplication_counts.txt",row.names=T,col.names=T,sep='\t',quote=F)

#Grouping 制作分组文件
#获取临床信息
GEO_clin<-data.frame("Sample"=colnames(Exp),"Type"=c(rep('UC',14),rep('control',16)))
Group<-GEO_clin
colnames(Group)<-c("Sample","Type")
Group<-Group[which(Group$Type %in% c('control','UC')),]

#排序
Group<-Group[order(Group$Type),]
table(Group$Type)
#control      UC 
#     16      14 

#保存分组文件
write.table(Group,'group.txt',col.names = TRUE,row.names = FALSE,sep = '\t',quote = FALSE)

#标准化
library(edgeR)
library(limma)
library(impute)

#data
#Exp<-read.table('deduplication_counts.txt',header = T,row.names = 1,sep = '\t')
exp<-Exp
exp<-exp[,Group$Sample]
dim(exp) #[1] 23777    30

#数据转换
#exp<-log2(exp+1)

#impute missing expression data
exp<-as.matrix(exp)
mat<-impute.knn(exp)
rt<-mat$data
dim(rt)
#[1] 23777    30

rt<-exp

boxplot(rt, col = "blue", xaxt = "n", outline = F)

#save plot
pdf(file="rawBox.pdf", width = 6, height = 6)
boxplot(rt, col = "blue", xaxt = "n", outline = F)
dev.off()


#normalize
rt2<-normalizeBetweenArrays(rt)
boxplot(rt2, col = "red", xaxt = "n", outline = F)
dim(rt2)
#[1] 23777    30

#save plot
pdf(file = "normalBox.pdf", width = 6, height = 6)
boxplot(rt2, col = "red", xaxt = "n", outline = F)
dev.off()

#save
write.table(rt2,"normalize.txt",row.names=T,quote=F,col.names=T,sep='\t')

