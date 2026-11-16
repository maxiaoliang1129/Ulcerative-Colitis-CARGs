

rm(list=ls())

#path
path<-'E:\\LZ-25007\\colitis\\06verify\\GSE87466\\00data'
setwd(path)

library(limma)
library(edgeR)
library(ggvenn)
library(patchwork)
library(data.table)

#data
data<-read.table("GSE87466_series_matrix.txt",header = T,row.names=1,sep = "\t")
Genesample<-fread("GPL13158-5065.txt",header = T,sep = "\t") #注释文件
dim(data);dim(Genesample)
#[1] 54715   108
#[1] 54715    16

head(Genesample)
colnames(Genesample)
Genesample<-Genesample[,c(1,11)]
colnames(Genesample)<-c('ID','Gene.Symbol')


#基因名
id<-Genesample$Gene.Symbol
for(i in 1:length(id)){
  a=strsplit(id[i],' /// ')[[1]][1]
  id[i]=a
}
Genesample$Gene.Symbol<-id
head(Genesample)

#去除gene Symbol为空的探针
Genesample<-na.omit(Genesample)
dim(Genesample) #[1] 41796     2

#----------------------------
data$ID<-rownames(data)
colnames(Genesample)[2]<-'Symbol'
mRNA_rt<-merge(data,Genesample,by='ID')
dim(mRNA_rt)
#[1] 41796   110

mRNA_rt<-mRNA_rt[,c('ID','Symbol',colnames(data))]

mRNA_rt<-mRNA_rt[,-ncol(mRNA_rt)]
mRNA_rt<-mRNA_rt[,-1]
colnames(mRNA_rt)[1]<-'gene'
mRNA_rt<-na.omit(mRNA_rt)
dim(mRNA_rt) #[1] 41796   109

#处理矩阵
merge_data<-mRNA_rt
gene<-data.frame(merge_data$gene)
exp<-merge_data[,2:ncol(merge_data)]
dim(gene)
dim(exp)
exp[1:4,1:3]

exp<-cbind(gene,exp)
colnames(exp)[1]<-"Symbol"
dim(exp)
exp[1:3,1:4]

rt<-exp
rt<-as.matrix(rt)
rownames(rt)<-rt[,1]
exp<-rt[,2:ncol(rt)]
dimnames<-list(rownames(exp),colnames(exp))
exp<-matrix(as.numeric(exp),nrow = nrow(exp), dimnames = dimnames)
dim(exp)  #[1] 54295   136

exp[1:3,1:3]

#数据转换
#exp<-log2(exp+1)

#Deduplication
counts<-exp
counts_delrep<-counts
counts_delrep<-avereps(counts_delrep)
Exp<-counts_delrep
dim(Exp)
#[1] 20277   108

write.table(Exp,"deduplication_counts.txt",row.names=T,col.names=T,sep='\t',quote=F)

#Grouping 制作分组文件
#获取临床信息
#获取临床信息
library(GEOquery)
GSE<-getGEO(filename="GSE87466_series_matrix.txt.gz",getGPL=F)
GSE_pheno<-phenoData(GSE)
GSE_pheno_df<-GSE_pheno@data
dim(GSE_pheno_df) #[1] 108  39

colnames(GSE_pheno_df)

cols<-c("geo_accession","disease:ch1")
GEO_clin<-GSE_pheno_df[,cols]
table(GEO_clin$`disease:ch1`)
#                 Normal Ulcerative Colitis (UC) 
#                     21                      87 

GEO_clin[which(GEO_clin$`disease:ch1` == 'Normal'),"disease:ch1"]<-'control'
GEO_clin[which(GEO_clin$`disease:ch1` == 'Ulcerative Colitis (UC)'),"disease:ch1"]<-'UC'
colnames(GEO_clin)<-c('Sample','Type')
table(GEO_clin$Type)
#control      UC 
#     21      87

Group<-GEO_clin
colnames(Group)<-c("Sample","Type")
Group<-Group[order(Group$Type),]

write.table(Group,"group.txt",row.names=F,quote=F,col.names=T,sep='\t')

#标准化
library(edgeR)
library(limma)
library(impute)

#data
#Exp<-read.table('deduplication_counts.txt',header = T,row.names = 1,sep = '\t')
exp<-Exp
exp<-exp[,Group$Sample]
dim(exp)

#数据转换
#exp<-log2(exp+1)

#impute missing expression data
exp<-as.matrix(exp)
mat<-impute.knn(exp)
rt<-mat$data
dim(rt)
#[1] 20277   108

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
#[1] 20277   108

#save plot
pdf(file = "normalBox.pdf", width = 6, height = 6)
boxplot(rt2, col = "red", xaxt = "n", outline = F)
dev.off()

#save
write.table(rt2,"normalize.txt",row.names=T,quote=F,col.names=T,sep='\t')

