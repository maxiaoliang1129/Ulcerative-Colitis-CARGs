

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\08ceRNA'
setwd(path)

#提取miRNA_mRNA文件
library(openxlsx)
mirtarbase<-read.xlsx("hsa_MTI.xlsx",sheet = 1)

'''
#筛选lncRNA相关的miRNA
gtf<-rtracklayer::import('gencode.v46.long_noncoding_RNAs.gtf')
gtf_df=as.data.frame(gtf)
mirtarbase_lncrna<-mirtarbase[which(mirtarbase[,4] %in% gtf_df$gene_name),]
#筛选TMPO-AS1:ENSG00000257167对应的miRNA
mirtarbase_lncrna_mirna<-mirtarbase_lncrna[which(mirtarbase_lncrna[,2] == 'TMPO-AS1'),]
'''
'''
lncRNA_miRNA<-read.table('miRNA.txt',header = T,sep = '\t')
TMPO.miRNA<-lncRNA_miRNA$miRNA.name
TMPO.miRNA<-unique(TMPO.miRNA)
length(TMPO.miRNA) #[1] 57
'''
'''
#从ENCORI数据库得到的结果
#hsa-miR-577
mirtarbase_TMPO<-mirtarbase[which(mirtarbase[,2] == 'hsa-miR-577'),]

#保留以构建网络
Net.result<-mirtarbase

#write.table(Net.result,'Net.result.txt',col.names = T,row.names = F,sep = '\t',quote = F)
'''


#提取mRNA
mRNA<-mirtarbase$Target.Gene
length(mRNA) #[1] 502652

mRNA<-unique(mRNA)
length(mRNA) #[1] 15064

#保存基因
#write.table(data.frame(mRNA),'mrna.gene.txt',col.names = T,row.names = F,sep = '\t',quote = F)

#输入基因集
PRGs<-read.table('venn.txt',header = F,sep = '\t')
PRGs<-PRGs$V1
length(PRGs) #[1] 2

#取交集
inter<-intersect(PRGs,mRNA)
length(inter) #[1] 2
inter
#[1] "CCL11" "CLDN1"


#构建miRNA-Gene网络
net1<-mirtarbase[which(mirtarbase$Target.Gene %in% inter),]
dim(net1)
#[1] 59  9

miRNA_mRNA<-net1[,c(2,4)]

#构建lncRNA-miRNA网络
#先提取miRNA
miRNA<-miRNA_mRNA$miRNA
length(miRNA) #[1] 59

#筛选lncRNA相关的miRNA
gtf<-rtracklayer::import('gencode.v46.long_noncoding_RNAs.gtf')
gtf_df=as.data.frame(gtf)
gtf_df<-gtf_df[which(gtf_df$gene_type == 'lncRNA'),]
dim(gtf_df)

mirtarbase_lncrna<-mirtarbase[which(mirtarbase[,4] %in% gtf_df$gene_name),]
dim(mirtarbase_lncrna)
#[1] 1851    9

#筛选miRNA对应的lncRNA
mirtarbase_lncrna_mirna<-mirtarbase_lncrna[which(mirtarbase_lncrna[,2] %in% miRNA),]
dim(mirtarbase_lncrna_mirna)
#[1] 96  9

net2<-mirtarbase_lncrna_mirna[,c(2,4)]
lncRNA_miRNA<-net2
dim(lncRNA_miRNA);dim(miRNA_mRNA)
#[1] 96  2
#[1] 59  2

#miRNA
a<-lncRNA_miRNA$miRNA
b<-miRNA_mRNA$miRNA
c<-intersect(a,b)
length(c)
#[1] 32

#去重
dim(lncRNA_miRNA) #[1] 52  2
lncRNA_miRNA<-unique(lncRNA_miRNA) #[1] 96  2

dim(miRNA_mRNA) #[1] 62  2
miRNA_mRNA<-unique(miRNA_mRNA) #[1] 59  2


network<-rbind(lncRNA_miRNA,miRNA_mRNA)

write.csv(network,'ceRNA.csv',row.names = F,quote = F)

'''
#########################数据222222222222
rm(list=ls())

#path
#path<-'D:\\LZ-24021\\lncRNA'
#setwd(path)

#提取miRNA_mRNA文件
library(openxlsx)
mirtarbase<-read.xlsx("hsa_MTI.xlsx",sheet = 1)

#筛选lncRNA相关的miRNA

'''
gtf<-rtracklayer::import('gencode.v46.long_noncoding_RNAs.gtf')
gtf_df=as.data.frame(gtf)
mirtarbase_lncrna<-mirtarbase[which(mirtarbase[,4] %in% gtf_df$gene_name),]
#筛选TMPO-AS1:ENSG00000257167对应的miRNA
mirtarbase_lncrna_mirna<-mirtarbase_lncrna[which(mirtarbase_lncrna[,2] == 'TMPO-AS1'),]
'''

lncRNA_miRNA<-read.table('miRNA.txt',header = T,sep = '\t')

#筛选lncRNA
lncRNA<-lncRNA_miRNA[which(lncRNA_miRNA$Experiments > 0 ),]
lncRNA<-lncRNA[which(lncRNA$Publications > 1 ),]
lncRNA<-lncRNA[which(lncRNA$Cell.Lines > 0 ),]
lncRNA<-lncRNA[which(lncRNA$Tissues > 0 ),]
dim(lncRNA) #[1] 12 7

TMPO.miRNA<-lncRNA$miRNA.name
TMPO.miRNA<-unique(TMPO.miRNA)
length(TMPO.miRNA) #[1] 12

TMPO.miRNA
'''
[1] "hsa-let-7a-5p"  "hsa-let-7b-5p"  "hsa-let-7d-5p"  "hsa-let-7e-5p"  "hsa-let-7g-5p" 
[6] "hsa-let-7i-5p"  "hsa-miR-26a-5p" "hsa-miR-26b-5p" "hsa-miR-30b-5p" "hsa-miR-30c-5p"
[11] "hsa-miR-96-5p"  "hsa-miR-98-5p" 
'''
#保存miRNA的结果
write.table(data.frame(TMPO.miRNA),'mirna.gene.txt',col.names = T,row.names = F,sep = '\t',quote = F)


#从ENCORI数据库得到的结果
mirtarbase_TMPO<-mirtarbase[which(mirtarbase[,2] %in% TMPO.miRNA),]

#保留以构建网络
Net.result<-mirtarbase_TMPO[,c(2,4)]

write.table(Net.result,'Net.result.txt',col.names = T,row.names = F,sep = '\t',quote = F)


#保留以构建网络
Net.result<-mirtarbase_TMPO[,c(2,4)]

write.table(Net.result,'Net.result.txt',col.names = T,row.names = F,sep = '\t',quote = F)


#提取mRNA
mRNA<-mirtarbase_TMPO$Target.Gene
length(mRNA) #[1] 9913

mRNA<-unique(mRNA)
length(mRNA) #[1] 4277

#保存基因
write.table(data.frame(mRNA),'mrna.gene.txt',col.names = T,row.names = F,sep = '\t',quote = F)

#输入基因集
PRGs<-read.table('diff_intergene.txt',header = F,sep = '\t')
PRGs<-PRGs$V1
length(PRGs) #[1] 15

#取交集
inter<-intersect(PRGs,mRNA)
length(inter) #[1] 4
inter
#[1] "CAT"     "NOS2"    "SLC27A2" "XDH"  
'''





