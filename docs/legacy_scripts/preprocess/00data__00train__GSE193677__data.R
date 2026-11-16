

rm(list=ls())

library(limma)
library(edgeR)
library(ggvenn)
library(patchwork)
library(data.table)

#path
path<-'E:\\LZ-25007\\colitis\\00data\\GSE193677'
setwd(path)

#data
data<-read.table("GSE193677_raw_counts_GRCh38.p13_NCBI.tsv",header = T,sep = "\t")
Genesample<-fread("Human.GRCh38.p13.annot.tsv",header = T,sep = "\t") #注释文件
dim(data);dim(Genesample)
#[1] 39376  2472
#[1] 39376    18

colnames(Genesample)
Genesample<-Genesample[,c(1,2)]
head(Genesample)

#基因名
id<-Genesample$Symbol
for(i in 1:length(id)){
  a=strsplit(id[i],' /// ')[[1]][1]
  id[i]=a
}
Genesample$Symbol<-id
head(Genesample)

#去除gene Symbol为空的探针
Genesample<-na.omit(Genesample)
dim(Genesample) #[1] 39376     2

#----------------------------
colnames(data)[1]<-'ID'
colnames(Genesample)<-c('ID','Symbol')
mRNA_rt<-merge(data,Genesample,by='ID')
dim(mRNA_rt)
#[1] 39376  2473

mRNA_rt<-mRNA_rt[,c('ID','Symbol',colnames(data)[2:ncol(data)])]

mRNA_rt<-mRNA_rt[,-1]
colnames(mRNA_rt)[1]<-'gene'
mRNA_rt<-na.omit(mRNA_rt)
dim(mRNA_rt) #[1] 39376  2472

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
dim(exp)  #[1] 39376  2471

exp[1:3,1:3]

#Deduplication
counts<-exp
counts_delrep<-counts
counts_delrep<-avereps(counts_delrep)
Exp<-counts_delrep
dim(Exp)
#[1] 39374  2471

write.table(Exp,"deduplication_counts.txt",row.names=T,col.names=T,sep='\t',quote=F)

#Grouping 制作分组文件
#获取临床信息
library(GEOquery)
GSE<-getGEO(filename="GSE193677_series_matrix.txt.gz",getGPL=F)
GSE_pheno<-phenoData(GSE)
GSE_pheno_df<-GSE_pheno@data
dim(GSE_pheno_df)

colnames(GSE_pheno_df)

#导出
write.csv(GSE_pheno_df,'GSE193677.csv')

cols<-c("geo_accession","ibd_disease:ch1")
GEO_clin<-GSE_pheno_df[,cols]
table(GEO_clin$`ibd_disease:ch1`)
#     CD Control      UC 
#   1157     461     872 

GEO_clin[which(GEO_clin$`ibd_disease:ch1` == 'Control'),"ibd_disease:ch1"]<-'control'
GEO_clin[which(GEO_clin$`ibd_disease:ch1` == 'UC'),"ibd_disease:ch1"]<-'UC'
colnames(GEO_clin)<-c('Sample','Type')
GEO_clin<-GEO_clin[which(GEO_clin$Type %in% c('control','UC')),]
table(GEO_clin$Type)
#control      UC 
#    461     872 

Group<-GEO_clin
colnames(Group)<-c("Sample","Type")
Group<-Group[order(Group$Type),]

write.table(Group,"group.txt",row.names=F,quote=F,col.names=T,sep='\t')

#标准化
library(edgeR)
library(limma)
library(impute)

#data
Exp<-read.table('deduplication_counts.txt',header = T,row.names = 1,sep = '\t')
exp<-Exp
#exp<-exp[,Group$Sample]

#表达矩阵和分组交集的样本数
a<-colnames(exp)
length(a)
#[1] 2471
b<-Group$Sample
length(b)
#[1] 1333
c<-intersect(a,b)
length(c)
#[1] 1320

#提取共同样本
exp<-exp[,c]
Group<-Group[which(Group$Sample %in% c),]
dim(Group) #[1] 1320    2
write.table(Group,"group2.txt",row.names=F,quote=F,col.names=T,sep='\t')

#数据转换
exp<-log2(exp+1)

#impute missing expression data
exp<-as.matrix(exp)
mat<-impute.knn(exp)
rt<-mat$data
dim(rt)
#[1] 39374  1320

boxplot(rt, col = "blue", xaxt = "n", outline = F)

#save plot
pdf(file="rawBox.pdf", width = 6, height = 6)
boxplot(rt, col = "blue", xaxt = "n", outline = F)
dev.off()


#normalize
rt2<-normalizeBetweenArrays(rt)
boxplot(rt2, col = "red", xaxt = "n", outline = F)

#save plot
pdf(file = "normalBox.pdf", width = 6, height = 6)
boxplot(rt2, col = "red", xaxt = "n", outline = F)
dev.off()

#save
write.table(rt2,"normalize.txt",row.names=T,quote=F,col.names=T,sep='\t')



################################统计唯一样本数##################################
library(data.table)
library(dplyr)

pheno <- fread("GSE193677.csv")

# 从 title 提取受试者 ID
pheno$patient_id <- sub(".*reGRID_([0-9]+)_Biopsy_.*", "\\1", pheno$title)

# 检查是否全部匹配
sum(!grepl("reGRID_[0-9]+_Biopsy_", pheno$title))  # 应为 0
#[1] 0

# 按疾病统计：活检数 vs 唯一患者数（含 CD，便于说明排除前结构）
pheno %>%
  group_by(`ibd_disease:ch1`) %>%
  summarise(
    n_biopsy = n(),
    n_patient = n_distinct(patient_id)
  )
'''
# A tibble: 3 × 3
  `ibd_disease:ch1` n_biopsy n_patient
  <chr>                <int>     <int>
1 CD                    1157       498
2 Control                461       243
3 UC                     872       421
'''

# 每个患者贡献几个活检的分布
pheno %>%
  count(patient_id, `ibd_disease:ch1`) %>%
  count(`ibd_disease:ch1`, n)

# 限制到最终纳入的 1,320 例（用 group2.txt）
group <- fread("group2.txt")
dat <- merge(pheno, group, by.x = "geo_accession", by.y = "Sample")

dat %>%
  group_by(Type) %>%
  summarise(
    n_biopsy = n(),
    n_patient = n_distinct(patient_id)
  )
'''
# A tibble: 2 × 3
  Type    n_biopsy n_patient
  <chr>      <int>     <int>
1 UC           865       419
2 control      455       241
'''
# UC 中 I/NonI 分布
dat %>%
  filter(Type == "UC") %>%
  count(`typere:ch1`)
'''
Key: <geo_accession>
   typere:ch1     n
       <char> <int>
1:          I   291
2:       NonI   574
'''
# UC 中 UC.I/UC.NonI 分布
dat %>%
  filter(Type == "UC") %>%
  count(`diseasetypere:ch1`)
'''
Key: <geo_accession>
   diseasetypere:ch1     n
              <char> <int>
1:              UC.I   291
2:           UC.NonI   574
'''


############################导出包含活检信息的分组信息##########################
library(data.table)
# dat 已由 pheno + group2 merge 得到

# 1. 同一患者是否只对应一个 Type（UC/control）——应为 0 行
dat[, .(n_type = uniqueN(Type)), by = patient_id][n_type > 1]

# 2. UC 患者内 I/NonI 是否一致——若 >0 行，说明亚型在活检间不稳定
dat[Type == "UC", .(n_sub = uniqueN(`diseasetypere:ch1`)), by = patient_id][n_sub > 1]

# 3. 限制到 1320 例后，每位患者贡献几个活检
dat[, .N, by = .(patient_id, Type)][, .N, by = .(Type, N)][order(Type, N)]
'''
      Type     N     N
    <char> <int> <int>
1:      UC     1   120
2:      UC     2   153
3:      UC     3   145
4:      UC     4     1
5: control     1    33
6: control     2   202
7: control     3     6
'''
dat[Type == "UC", .(n_patient = uniqueN(patient_id)), by = `diseasetypere:ch1`]
'''
   diseasetypere:ch1 n_patient
              <char>     <int>
1:           UC.NonI       375
2:              UC.I       229
'''

fwrite(dat[, .(Sample = geo_accession, Type, Patient = patient_id)],
       "sample_patient_group.txt", sep = "\t")

######################敏感性分析用：每位患者一个活检（两列）####################
set.seed(20240101)
one <- dat[, .SD[sample(.N, 1)], by = patient_id]
fwrite(one[, .(Sample = geo_accession, Type)],
       "group_one_per_patient.txt", sep = "\t")

dim(one)                 # 应为 660 行
#[1] 660  83
table(one$Type)          # UC 419 / control 241
#control      UC 
#  241     419 


#####这个是首次活检，这里选择使用随机活检就好了
dat[, biopsy_num := as.integer(sub(".*_Biopsy_([0-9]+).*", "\\1", title))]
one <- dat[order(patient_id, biopsy_num), .SD[1], by = patient_id]



