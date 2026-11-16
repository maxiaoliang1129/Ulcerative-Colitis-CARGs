

rm(list=ls())

path<-'D:\\LZ-25234\\UC\\12Subtype\\01ConsensusClusterPlus'
setwd(path)

library(ConsensusClusterPlus)

#data
data<-read.table('normalize.txt',header = T,sep='\t')
group<-read.table('group.txt',header = T,sep = '\t')
group<-group[which(group$Type == 'UC'),]

Exp<-data[,group$Sample]
dim(Exp)
#[1] 39374  865

inter_gene<-read.table('venn.txt',header=F)
Exp<-Exp[inter_gene$V1,]
dim(Exp) #[1]  2 865
rownames(Exp)<-toupper(rownames(Exp))

###聚类
df<-Exp

#这里没有进行中位数归一化
df<-as.matrix(df)
#进行一致性聚类
ccres<-ConsensusClusterPlus(df,maxK=9,reps=100,pItem=0.8,pFeature=1,
                            tmyPal = c("white","#C75D30"),
                            title='ConsensusCluster/',
                            clusterAlg="km",
                            distance="euclidean",
                            seed=123456,
                            plot="pdf")
iclres <- calcICL(ccres,title="ics of ssgsea res")

#根据PAC确定最佳聚类数
##PAC = Proportion of ambiguous clustering 模糊聚类比例
Kvec=2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 

names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = ccres[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK #[1] 2

#提取结果
sample_subtypes <- ccres[[2]][["consensusClass"]]
table(sample_subtypes)
#sample_subtypes
#sample_subtypes
#1  2 
#474 391

group<-sample_subtypes
NMF_result <- as.data.frame(group)
head(NMF_result,2)

write.csv(NMF_result,"NMF_result.csv")


#PCA plot2
table(group)
#group
#  1   2 
# 474 391

group<-as.factor(group)
pca_result<-prcomp(t(df),scale=T)

library(ggbiplot)
ggbiplot(pca_result,
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,          # 横纵比例
         groups = group,     # 添加分组信息，为分组文件的第一列
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+ylim(c(-3,3))

#plot
pdf(file="pca.pdf",width=4,height=4)
ggbiplot(pca_result,
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,          # 横纵比例
         groups = group,     # 添加分组信息，为分组文件的第一列
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+ylim(c(-3,3))
dev.off()

'lightblue''steelblue3'
'lightpink''cadetblue1'



################################补充分析1#######################################
library(cluster)
library(mclust)

# --- silhouette（与聚类一致用 euclidean）---
dist_mat <- dist(t(df), method = "euclidean")
sil <- silhouette(as.numeric(ccres[[2]]$consensusClass), dist_mat)
avg_sil <- summary(sil)$avg.width
cat("Average silhouette width for k=2:", avg_sil, "\n")
#Average silhouette width for k=2: 0.3622578

sil_df <- as.data.frame(sil[, 1:3])
colnames(sil_df) <- c("cluster", "neighbor", "sil_width")
write.csv(sil_df, "silhouette_k2.csv", row.names = FALSE)

pdf("silhouette_k2.pdf", width = 6, height = 4)
plot(sil, col = c("lightpink", "cadetblue1"), main = "Silhouette k=2")
dev.off()

# --- calcICL 的 clusterConsensus 作为稳定性指标 ---
icl <- calcICL(ccres, title = "ICL", plot = "pdf")
write.csv(icl$clusterConsensus, "cluster_consensus.csv", row.names = FALSE)


#################################补充分析2######################################
# --- 留两个基因聚类 ---
# full_exp 是全部 UC 样本的完整表达矩阵
full_exp <- data[, group$Sample]
rownames(full_exp) <- toupper(rownames(full_exp))   # 确保与 CCL11/CLDN1 写法一致

keep_genes <- setdiff(rownames(full_exp), c("CCL11", "CLDN1"))
lofo_exp <- full_exp[keep_genes, ]

# 关键：强制转成 matrix
lofo_exp <- as.matrix(lofo_exp)
storage.mode(lofo_exp) <- "numeric"   # 防止有字符列导致后续报错

lofo_cc <- ConsensusClusterPlus(
  lofo_exp,
  maxK = 2, reps = 100, pItem = 0.8, pFeature = 1,
  tmyPal = c("white", "#C75D30"),
  title = 'ConsensusCluster_LOFO/',
  clusterAlg = "pam",        # 改用 pam 算法
  distance = "euclidean",
  seed = 123456,
  plot = "pdf"
)
lofo_class <- lofo_cc[[2]]$consensusClass
orig_class <- ccres[[2]]$consensusClass

ari <- adjustedRandIndex(orig_class, lofo_class)
cat("ARI (original vs LOFO):", ari, "\n")
write.csv(data.frame(ARI = ari), "ARI_LOFO.csv", row.names = FALSE)

# 置换检验给 ARI 一个 p 值（可选但推荐）
set.seed(123)
perm_ari <- replicate(200, adjustedRandIndex(sample(orig_class), lofo_class))
p_ari <- mean(perm_ari >= ari)
cat("Permutation p for ARI:", p_ari, "\n")
write.csv(data.frame(ARI = ari, p_perm = p_ari), "ARI_LOFO_perm.csv", row.names = FALSE)

#######输出
sample_subtypes <- ccres[[2]][["consensusClass"]]
# 改为 Group 标签
group_label <- ifelse(sample_subtypes == 1, "Group1", "Group2")
NMF_result <- data.frame(id = names(group_label), Group = group_label)
write.csv(NMF_result, "NMF_result.csv", row.names = FALSE)





rm(list=ls())

path <- 'D:\\LZ-25234\\UC\\12Subtype\\01ConsensusClusterPlus'
setwd(path)

library(ConsensusClusterPlus)
library(cluster)
library(mclust)

# ============================================================
# 1. 数据准备
# ============================================================
data <- read.table('normalize.txt', header = TRUE, sep = '\t')

# 兼容处理：若第一列是基因名（字符），设为行名并删除该列
if (is.character(data[[1]]) || is.factor(data[[1]])) {
  rownames(data) <- toupper(as.character(data[[1]]))
  data <- data[, -1, drop = FALSE]
} else {
  rownames(data) <- toupper(rownames(data))
}

group_df <- read.table('group.txt', header = TRUE, sep = '\t')
group_df <- group_df[which(group_df$Type == 'UC'), ]

# 关键：用独立变量保存 UC 样本 ID，避免后续 group 被覆盖
uc_samples <- group_df$Sample

# UC 样本表达矩阵（完整转录组）
full_exp <- as.matrix(data[, uc_samples])
storage.mode(full_exp) <- "numeric"
dim(full_exp)   # 应为 39374 x 865

# 两个基因的表达矩阵（用于原始聚类）
inter_gene <- read.table('venn.txt', header = FALSE)
Exp2 <- full_exp[inter_gene$V1, , drop = FALSE]
dim(Exp2)       # 应为 2 x 865

# 剔除 CCL11/CLDN1 后的完整转录组（用于 PCA、silhouette、LOFO）
keep_genes <- setdiff(rownames(full_exp), c("CCL11", "CLDN1"))
lofo_exp <- full_exp[keep_genes, , drop = FALSE]
dim(lofo_exp)   # 应为 39372 x 865

# ============================================================
# 2. 原始共识聚类（基于 CCL11 和 CLDN1）
# ============================================================
df <- as.matrix(Exp2)
storage.mode(df) <- "numeric"

ccres <- ConsensusClusterPlus(
  df, maxK = 9, reps = 100, pItem = 0.8, pFeature = 1,
  tmyPal = c("white", "#C75D30"),
  title = 'ConsensusCluster/',
  clusterAlg = "km",
  distance = "euclidean",
  seed = 123456,
  plot = "pdf"
)

# PAC 确定最佳 k
Kvec <- 2:9
x1 <- 0.1; x2 <- 0.9
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste("K=", Kvec, sep = "")
for (i in Kvec) {
  M <- ccres[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
}
optK <- Kvec[which.min(PAC)]
cat("Optimal K by PAC:", optK, "\n")
#Optimal K by PAC: 2 
write.csv(data.frame(K = Kvec, PAC = PAC), "PAC.csv", row.names = FALSE)

# 提取原始分组，统一为 Group1/Group2
sample_subtypes <- ccres[[2]][["consensusClass"]]
table(sample_subtypes)
'''
sample_subtypes
  1   2 
474 391 
'''

group_label <- ifelse(sample_subtypes == 1, "Group1", "Group2")
names(group_label) <- names(sample_subtypes)
NMF_result <- data.frame(id = names(group_label), Group = group_label, stringsAsFactors = FALSE)
write.csv(NMF_result, "NMF_result.csv", row.names = FALSE)

# ============================================================
# 3. 稳定性评估
# ============================================================

# 3.1 silhouette —— 基于完整转录组（剔除 CCL11/CLDN1，避免循环论证）
dist_full <- dist(t(lofo_exp), method = "euclidean")
sil <- silhouette(as.numeric(sample_subtypes), dist_full)
avg_sil <- summary(sil)$avg.width
cat("Average silhouette width (full transcriptome, excluding CCL11/CLDN1):", avg_sil, "\n")
#Average silhouette width (full transcriptome, excluding CCL11/CLDN1): 0.04181778 
sil_df <- as.data.frame(sil[, 1:3])
colnames(sil_df) <- c("cluster", "neighbor", "sil_width")
write.csv(sil_df, "silhouette_full_transcriptome.csv", row.names = FALSE)

pdf("silhouette_full_transcriptome.pdf", width = 6, height = 4)
plot(sil, col = c("lightpink", "cadetblue1"),
     main = "Silhouette (full transcriptome, CCL11/CLDN1 excluded)")
dev.off()

# 3.2 calcICL 共识值
icl <- calcICL(ccres, title = "ICL", plot = "pdf")
write.csv(icl$clusterConsensus, "cluster_consensus.csv", row.names = FALSE)

# 3.3 PCA —— 基于完整转录组（剔除 CCL11/CLDN1），避免循环论证
pca_result <- prcomp(t(lofo_exp), scale = TRUE)
library(ggbiplot)
pca_group <- factor(group_label, levels = c("Group1", "Group2"))

p_pca <- ggbiplot(pca_result,
                  var.axes = FALSE,
                  obs.scale = 1,
                  groups = pca_group,
                  ellipse = TRUE,
                  circle = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank()) +
  labs(color = "CCL11/CLDN1 group") +
  ggtitle("PCA on full transcriptome (CCL11/CLDN1 excluded)")

pdf("pca_full_transcriptome.pdf", width = 4, height = 4)
print(p_pca)
dev.off()

# ============================================================
# 4. Leave-two-genes-out 验证 + ARI（方案二：hclust + cutree）
# ============================================================
# 用剔除 CCL11/CLDN1 后的完整转录组，直接进行层次聚类，切成 2 类
dist_lofo <- dist(t(lofo_exp), method = "euclidean")
hc_lofo <- hclust(dist_lofo, method = "average")
lofo_class_direct <- cutree(hc_lofo, k = 2)

# 对齐样本顺序后计算 ARI
orig_class <- sample_subtypes
lofo_class_direct <- lofo_class_direct[names(orig_class)]

ari_direct <- adjustedRandIndex(orig_class, lofo_class_direct)
cat("ARI (direct hclust, LOFO):", ari_direct, "\n")
#ARI (direct hclust, LOFO): -0.000404853 

# 置换检验
set.seed(123)
perm_ari <- replicate(200, adjustedRandIndex(sample(orig_class), lofo_class_direct))
p_ari <- mean(perm_ari >= ari_direct)
cat("Permutation p for ARI:", p_ari, "\n")
#Permutation p for ARI: 1

write.csv(data.frame(ARI = ari_direct, p_perm = p_ari),
          "ARI_LOFO_direct.csv", row.names = FALSE)

# 可选：把 LOFO 分组也保存下来，方便后续敏感性分析
lofo_group_label <- ifelse(lofo_class_direct == 1, "LOFO_Group1", "LOFO_Group2")
lofo_df <- data.frame(id = names(lofo_group_label),
                      LOFO_Group = lofo_group_label,
                      stringsAsFactors = FALSE)
write.csv(lofo_df, "LOFO_group_assignment.csv", row.names = FALSE)


# ============================================================
# 5. 可选：bootstrap Jaccard（若已安装 fpc 包；大数据较慢，可跳过）
# ============================================================
if (requireNamespace("fpc", quietly = TRUE)) {
  library(fpc)
  set.seed(123)
  # 为加速，仅用方差最大的前 2000 个基因做 bootstrap
  gene_var <- apply(lofo_exp, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, length(gene_var))]
  boot_mat <- t(lofo_exp[top_genes, , drop = FALSE])
  
  boot_res <- clusterboot(boot_mat,
                          clustermethod = kmeansCBI,
                          k = 2, B = 100)
  cat("Bootstrap Jaccard (LOFO, top 2000 variable genes):\n")
  print(boot_res$bootmean)
  write.csv(data.frame(cluster = 1:2, Jaccard = boot_res$bootmean),
            "bootstrap_jaccard.csv", row.names = FALSE)
} else {
  message("fpc package not installed; skipping bootstrap Jaccard.")
}
#fpc package not installed; skipping bootstrap Jaccard.



rm(list=ls())

library(ggplot2)
library(ggpubr)

#check points gene表达与绘图
#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\02boxplot'
setwd(path)

#特征基因表达文件
rt<-read.table('normalize.txt',sep = '\t',header = T,row.names = 1)
group<-read.table('group.txt',header = T,sep = '\t')
group<-group[which(group$Type == 'UC'),]
rt<-rt[,group$Sample]
dim(rt)

check_gene<-read.table('venn.txt',sep = '\t')

rows<-na.omit(match(check_gene$V1,rownames(rt)))
rt_imm_checkp<-rt[rows,]

rownames(rt_imm_checkp)<-toupper(rownames(rt_imm_checkp))

write.table(rt_imm_checkp, file="hub_exp.txt",
            sep="\t",quote=F, col.names=T, row.names=T)

#boxplot
rt<-read.table("hub_exp.txt",sep="\t",header=T,row.names=1,check.names=F)
Type<-read.csv("NMF_result.csv",sep = ",",header = T)

Type<-Type[order(Type$group),]
rt<-rt[,Type$X]

cols<-na.omit(match(rownames(Type), colnames(rt)))
rt=t(rt)
rownames(Type)<-Type$X
rownames(rt)==rownames(Type)

data = data.frame()
rows = na.omit(match(rownames(rt), rownames(Type)))
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[rows, 2])))
  a[ , 1] = as.numeric(a[ , 1])
  test <- kruskal.test(expression ~ Subtype, data = a)
  pValue = test$p.value
  if(!is.na(pValue) & pValue < 0.05){
    data = rbind(data, a)
    print(pValue)
  }
}

data[which(data$Subtype == 1),'Subtype']<-'Subtype1'
data[which(data$Subtype == 2),'Subtype']<-'Subtype2'
#data[which(data$Subtype == 'cluster3'),'Subtype']<-'cluster3'
data<-data[order(data$Subtype, decreasing = T), ]

write.table(data, file="P0.05_data.txt",sep="\t",row.names=F,quote=F)

'''
#plot_0.01
data=read.table("P0.05_data.txt",sep="\t",header=T,check.names=F)  
#data$expression = log2(data$expression + 0.01)
data$Subtype = factor(data$Subtype, levels = c("Subtype1", "Subtype2"))
p = ggboxplot(data, x="gene", y="expression", fill = "Subtype",
              ylab = "Gene expression",
              xlab = "", palette = c("red", "blue"))+
  theme(text = element_text(size = 21))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Subtype),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif")
p
ggsave(filename = "boxplot1.tiff", plot = p, dpi = 400, width = 17, height = 8)
ggsave(filename = "boxplot1.pdf", plot = p, dpi = 500, width = 17, height = 8)
'''

#plot_all
data = data.frame()
rows = na.omit(match(rownames(rt), rownames(Type)))
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[rows, 2])))
  a[ , 1] = as.numeric(a[ , 1])
  test <- kruskal.test(expression ~ Subtype, data = a)
  pValue = test$p.value
  if(!is.na(pValue)){
    data = rbind(data, a)
    print(pValue)
  }
}

data[which(data$Subtype == 1),'Subtype']<-'Subtype1'
data[which(data$Subtype == 2),'Subtype']<-'Subtype2'
#data[which(data$Subtype == 'cluster3'),'Subtype']<-'cluster3'
data<-data[order(data$Subtype, decreasing = T), ]

write.table(data, file="Pall_data.txt",sep="\t",row.names=F,quote=F)

#plot
data=read.table("Pall_data.txt",sep="\t",header=T,check.names=F)  
#data$expression = log2(data$expression + 0.01)
data$Subtype = factor(data$Subtype, levels = c("Subtype1", "Subtype2"))
p = ggboxplot(data, x="gene", y="expression", fill = "Subtype",
              ylab = "Gene expression",
              xlab = "", palette = c('lightpink','cadetblue1'))+
  theme(text = element_text(size = 12))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Subtype),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif")
p

ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 3, height = 4)





rm(list=ls())

#TCGA免疫微环境分析和免疫治疗反应性预测

library(GSVA)
library(limma)
library(GSEABase)

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\03immue\\01data'
setwd(path)

#data
rt<-read.table('normalize.txt',sep = '\t',header = T,row.names = 1)
Risk_sample<-read.csv('NMF_result.csv',header = TRUE,sep = ',',check.names=F)
colnames(Risk_sample)[1]<-'id'
dim(rt);dim(Risk_sample)
#[1] 39374  1320
#[1] 865   2

rt<-rt[,Risk_sample$id]
dim(rt) #[1] 39374   865

'''
for(i in 1:length(Risk_sampleid)){
  Risk_sampleid[i]<-gsub('-','.',Risk_sampleid[i])
}
rt<-rt[,Risk_sampleid]
dim(rt)
'''

exp<-rt
dimnames=list(rownames(exp), colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0, ]
dim(mat)
#[1] 39374   865

gmtFile="immune.gmt"

'''
x<-readLines(gmtFile)
res<-strsplit(x, "\t")
names(res)<-vapply(res, function(y) y[1],character(1))
res<-lapply(res,"[",-c(1:2))
'''

geneSet=getGmt(gmtFile, geneIdType = SymbolIdentifier())

ssgseaScore=gsva(ssgseaParam(exprData=mat,geneSets = geneSet))

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))
}

ssgseaOut = normalize(ssgseaScore)
ssgseaOut = rbind(id=colnames(ssgseaOut),ssgseaOut)

write.table(ssgseaOut, file = "ssgseaOut.txt", sep = "\t", quote = F, col.names = F)




rm(list=ls())

#estimate计算分值

library(limma)
library(estimate)

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\03immue\\02estimate'
setwd(path)

#data
rt<-read.table('normalize.txt',sep = '\t',header = T,row.names = 1)
Risk_sample<-read.csv('NMF_result.csv',header = TRUE,sep = ',',check.names=F)
colnames(Risk_sample)[1]<-'id'
dim(rt);dim(Risk_sample)
#[1] 39374  1320
#[1] 865   2

rt<-rt[,Risk_sample$id]

'''
Risk_sampleid<-Risk_sample$id
for(i in 1:length(Risk_sampleid)){
  Risk_sampleid[i]<-gsub('-','.',Risk_sampleid[i])
}
rt<-rt[,Risk_sampleid]
dim(rt)
'''

exp = rt
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
dim(data)
#[1] 39374   865

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

 

rm(list=ls())

#estimate heatmap
library(ggsci)
library(pheatmap)

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\03immue\\03estimateHM'
setwd(path)

#data
rt<-read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)
Risk_sample<-read.csv('NMF_result.csv',header = TRUE,sep = ',',check.names=F)
colnames(Risk_sample)[1]<-'id'
clust<-Risk_sample
clust$cluster<-'cluster'
clust[which(clust$group == 1),'cluster']<-'cluster1'
clust[which(clust$group == 2),'cluster']<-'cluster2'
#clust[which(clust$group == 3),'cluster']<-'cluster3'
Type<-clust
dim(clust);head(clust)

write.table(clust,'cluster.txt',sep = '\t',quote = F,row.names = F)

score=read.table("scores.txt",sep="\t",check.names=F, header=T)
colnames(score)[1]<-'id'
score$id<-gsub('-','.',score$id)
cluster = merge.data.frame(Type, score, by = "id")

#rt = rt[,cluster[,'id']]
cluster<-cluster[order(cluster$cluster),]
rt<-rt[,cluster$id]

unique(cluster$V1 == colnames(rt))  #check
rownames(cluster) = cluster$id
cluster = cluster[ , -c(1, 2)]
colnames(cluster)[1] = "Subtype"
cluster = cluster[order(cluster$Subtype),]


#plot
cluster[which(cluster$Subtype == 'cluster1'),'Subtype']<-'Subtype1'
cluster[which(cluster$Subtype == 'cluster2'),'Subtype']<-'Subtype2'
#cluster[which(cluster$Subtype == 'cluster3'),'Subtype']<-'Subtype3'
cluster<-cluster[,-5]

ann_color <- list(Subtype = c(Subtype1 ='lightpink', Subtype2 = 'cadetblue1'))

pheatmap(rt, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)

pdf("estimateHM.pdf", height = 4.5, width = 6.5,onefile = T)
pheatmap(rt, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\03immue\\03estimateHM'
setwd(path)

rt=read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)   #¶ÁÈ¡ÎÄ¼þ
Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=T)

rt=t(rt[,row.names(Type)])
unique(rownames(rt) == rownames(Type))
colnames(rt)

num_immu_cell=c(1,4,6,9,11,13,14,16,17,19,22,23,24,25,26,27)
num_immu_func = setdiff(1:length(colnames(rt)), num_immu_cell)

#num_immu_func
data = data.frame()
for(i in colnames(rt)[num_immu_func]){
  a = as.data.frame(cbind(expression=rt[,i], gene=i, Subtype=as.vector(Type[,1])))
  a[,1] = as.numeric(a[,1])
  wilcoxon <- t.test(expression ~ Subtype, data = a)
  pValue = wilcoxon$p.value
#  if(!is.na(pValue) & pValue < 0.05){
#    data = rbind(data, a)
  print(pValue)
#  }
  data=rbind(data, cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[,1])))
}

#save
write.table(data,file="data_func.txt",sep="\t",row.names=F,quote=F)


#num_immu_cell
data = data.frame()
for(i in colnames(rt)[num_immu_cell]){
  a = as.data.frame(cbind(expression=rt[,i], gene=i, Subtype=as.vector(Type[,1])))
  a[,1] = as.numeric(a[,1])
  wilcoxon <- t.test(expression ~ Subtype, data = a)
  pValue = wilcoxon$p.value
#  if(!is.na(pValue) & pValue < 0.05){
#    data = rbind(data, a)
  print(pValue)
#  }
  data=rbind(data, cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[,1])))
}

#save
write.table(data,file="data_cell.txt",sep="\t",row.names=F,quote=F)

#画图
library(ggplot2)
library(ggpubr)

#plot num_immu_cell
data=read.table("data_cell.txt",sep="\t",header=T,check.names=F)  
data$Subtype = factor(data$Subtype,levels=c(1,2))
colnames(data)[3]<-'Group'

p = ggboxplot(data, x="gene", y="expression", fill = "Group",
              ylab = "Fraction",
              xlab = "", palette = c("lightpink","cadetblue1"))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 5)+
  theme(text = element_text(size = 12))
p

ggsave(filename = "immu_cell_boxplot.pdf", plot = p, width = 5, height = 5)


#plot num_immu_func
data=read.table("data_func.txt",sep="\t",header=T,check.names=F)  
data$Subtype = factor(data$Subtype,levels=c(1,2))
colnames(data)[3]<-'Group'

p = ggboxplot(data, x="gene", y="expression", fill = "Group",
              ylab = "Fraction",
              xlab = "", palette = c("lightpink","cadetblue1"))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 5)+
  theme(text = element_text(size = 12))
p

ggsave(filename = "immu_func_boxplot.pdf", plot = p, width = 5, height = 6)

#合并绘图
data1=read.table("data_cell.txt",sep="\t",header=T,check.names=F)  
data2=read.table("data_func.txt",sep="\t",header=T,check.names=F) 
data<-rbind(data1,data2)

data$Subtype = factor(data$Subtype,levels=c(1,2))
colnames(data)[3]<-'Group'

p = ggboxplot(data, x="gene", y="expression", fill = "Group",
              ylab = "Fraction",
              xlab = "", palette = c("lightpink","cadetblue1"))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 5)+
  theme(text = element_text(size = 12))
p

ggsave(filename = "all_boxplot.pdf", plot = p, width = 8,
       height = 4.5)

###############################################################################
## ssGSEA / ESTIMATE 评分 GEE：调整 Activity + Endoscopy + Inflammation + BiopsySite
## 输入：ssgseaOut.txt、scores.txt、NMF_result.csv、clinical_metadata.txt
## 输出：GEE_ssgsea_adjusted.csv、GEE_estimate_adjusted.csv
###############################################################################

rm(list = ls())
setwd('D:\\LZ-25234\\UC\\12Subtype\\03immue\\03estimateHM')

if (!requireNamespace("geepack", quietly = TRUE)) install.packages("geepack")
library(geepack)

options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 1. 读分组
## ---------------------------------------------------------------------------
grp <- read.csv("NMF_result.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(grp)[1:2] <- c("id", "Group")
grp$id <- as.character(grp$id)
grp$Group <- as.character(grp$Group)
grp$Group[grp$Group == "1"] <- "Group1"
grp$Group[grp$Group == "2"] <- "Group2"
grp$Group <- factor(grp$Group, levels = c("Group1", "Group2"))
if (any(is.na(grp$Group))) stop("NMF_result.csv 的 Group 只能包含 1/2 或 Group1/Group2")

## ---------------------------------------------------------------------------
## 2. 读临床表（统一因子化，避免 character 报错）
## ---------------------------------------------------------------------------
clin <- read.table("clinical_metadata.txt",
                   sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
clin$id <- as.character(clin$id)
clin$Activity     <- factor(clin$Activity, levels = c("Inactive", "Active"))
clin$Inflammation <- factor(clin$Inflammation, levels = c("NonInflamed", "Inflamed"))
clin$BiopsySite   <- factor(clin$BiopsySite)
clin$Patient_ID   <- factor(clin$Patient_ID)
clin$Endoscopy    <- factor(clin$Endoscopy,
                            levels = sort(unique(clin$Endoscopy[!is.na(clin$Endoscopy)])))

## ---------------------------------------------------------------------------
## 3. 读 ssGSEA 评分
## 注意：ssgseaOut.txt 实际是 行=免疫条目，列=样本，需要转置
## ---------------------------------------------------------------------------
ssg_raw <- read.table("ssgseaOut.txt",
                      sep = "\t", header = TRUE, row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)

## 自动判断：如果列名是 GSM 样本名，说明需要转置
if (all(grepl("^GSM", colnames(ssg_raw)))) {
  ssg <- data.frame(t(ssg_raw), check.names = FALSE)
  ssg$id <- rownames(ssg)
  rownames(ssg) <- NULL
} else {
  ssg <- ssg_raw
  ssg$id <- rownames(ssg)
  rownames(ssg) <- NULL
}

score_cols_ssg <- setdiff(colnames(ssg), "id")
cat("ssGSEA 评分条目数：", length(score_cols_ssg), "\n")
#ssGSEA 评分条目数： 29 
cat("ssGSEA 样本数：", nrow(ssg), "\n")
#ssGSEA 样本数： 865 

## 检查 id 是否匹配，便于排查
cat("ssGSEA id 示例：", paste(head(ssg$id), collapse = ", "), "\n")
cat("分组 id 示例：", paste(head(grp$id), collapse = ", "), "\n")
cat("ssGSEA 与分组 id 交集数：", length(intersect(ssg$id, grp$id)), "\n")
#ssGSEA 与分组 id 交集数： 865

## 检查全 0 / 全 NA 列
bad_ssg <- score_cols_ssg[
  apply(ssg[, score_cols_ssg, drop = FALSE], 2,
        function(x) all(is.na(x)) || all(x == 0, na.rm = TRUE))
]
if (length(bad_ssg) > 0) {
  warning("以下 ssGSEA 条目全为 0/NA，将剔除：", paste(bad_ssg, collapse = ", "))
  score_cols_ssg <- setdiff(score_cols_ssg, bad_ssg)
}


## ---------------------------------------------------------------------------
## 4. 读 ESTIMATE 分数（行=样本，列=StromalScore/ImmuneScore/ESTIMATEScore）
## ---------------------------------------------------------------------------
est <- read.table("scores.txt",
                  sep = "\t", header = TRUE, row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
est$id <- rownames(est)

## 你原脚本里有 gsub('-', '.', score$id)，说明 scores.txt 的 id 可能用点号。
## 这里统一成 GSM 号：如果 id 是 "GSM5976499" 就不用改；
## 如果是 "GSM.5976499" 这种，请按实际情况调整下面这行。
# est$id <- gsub("\\.", "", est$id)  # 仅当 id 被点号拆分时使用，按需取消注释

score_cols_est <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
                            colnames(est))
if (length(score_cols_est) == 0) {
  warning("scores.txt 未找到 StromalScore/ImmuneScore/ESTIMATEScore，将使用除 id 外所有数值列")
  score_cols_est <- setdiff(colnames(est), "id")
}

## ---------------------------------------------------------------------------
## 5. 合并 + 完整案例集
## ---------------------------------------------------------------------------
dat_ssg <- merge(ssg,  grp,  by = "id")
dat_ssg <- merge(dat_ssg, clin, by = "id")

dat_est <- merge(est,  grp,  by = "id")
dat_est <- merge(dat_est, clin, by = "id")

gee_vars <- c("Group", "Activity", "Endoscopy", "Inflammation", "BiopsySite", "Patient_ID")

prep_cc <- function(dat, score_cols) {
  dat_cc <- dat[complete.cases(dat[, c(score_cols, gee_vars), drop = FALSE]), ]
  dat_cc$Group        <- droplevels(factor(dat_cc$Group, levels = c("Group1", "Group2")))
  dat_cc$Patient_ID   <- droplevels(factor(dat_cc$Patient_ID))
  dat_cc$Activity     <- droplevels(factor(dat_cc$Activity, levels = c("Inactive", "Active")))
  dat_cc$Endoscopy    <- droplevels(factor(dat_cc$Endoscopy))
  dat_cc$Inflammation <- droplevels(factor(dat_cc$Inflammation, levels = c("NonInflamed", "Inflamed")))
  dat_cc$BiopsySite   <- droplevels(factor(dat_cc$BiopsySite))
  dat_cc
}

dat_ssg_cc <- prep_cc(dat_ssg, score_cols_ssg)
dat_est_cc <- prep_cc(dat_est, score_cols_est)

cat("ssGSEA 完整案例样本数：", nrow(dat_ssg_cc),
    "，患者数：", length(unique(dat_ssg_cc$Patient_ID)), "\n")
#ssGSEA 完整案例样本数： 860 ，患者数： 416 

cat("ESTIMATE 完整案例样本数：", nrow(dat_est_cc),
    "，患者数：", length(unique(dat_est_cc$Patient_ID)), "\n")
#ESTIMATE 完整案例样本数： 860 ，患者数： 416 

if (nlevels(dat_ssg_cc$Group) < 2) stop("ssGSEA 完整案例中 Group 只有一个水平")
if (nlevels(dat_est_cc$Group) < 2) stop("ESTIMATE 完整案例中 Group 只有一个水平")

## ---------------------------------------------------------------------------
## 6. GEE 函数：未调整 vs 调整
## ---------------------------------------------------------------------------
run_gee <- function(dat_cc, score_cols, label) {
  out <- data.frame()
  for (s in score_cols) {
    f0 <- as.formula(paste0("`", s, "` ~ Group"))
    f1 <- as.formula(paste0("`", s, "` ~ Group + Activity + Endoscopy + Inflammation + BiopsySite"))
    
    m0 <- tryCatch(geeglm(f0, data = dat_cc, id = Patient_ID,
                          family = gaussian, corstr = "exchangeable"),
                   error = function(e) NULL)
    m1 <- tryCatch(geeglm(f1, data = dat_cc, id = Patient_ID,
                          family = gaussian, corstr = "exchangeable"),
                   error = function(e) NULL)
    
    if (is.null(m0) || is.null(m1)) {
      out <- rbind(out, data.frame(
        Score = s, Beta_unadj = NA, P_unadj = NA,
        Beta_adj = NA, P_adj = NA, Attenuation = NA,
        N = nrow(dat_cc), Note = "模型拟合失败"
      ))
      next
    }
    
    coef0 <- summary(m0)$coefficients
    coef1 <- summary(m1)$coefficients
    if (!"GroupGroup2" %in% rownames(coef0) || !"GroupGroup2" %in% rownames(coef1)) {
      out <- rbind(out, data.frame(
        Score = s, Beta_unadj = NA, P_unadj = NA,
        Beta_adj = NA, P_adj = NA, Attenuation = NA,
        N = nrow(dat_cc), Note = "GroupGroup2 系数缺失"
      ))
      next
    }
    
    b0 <- coef0["GroupGroup2", "Estimate"]
    p0 <- coef0["GroupGroup2", "Pr(>|W|)"]
    b1 <- coef1["GroupGroup2", "Estimate"]
    p1 <- coef1["GroupGroup2", "Pr(>|W|)"]
    atten <- ifelse(abs(b0) > 1e-4, 1 - abs(b1) / abs(b0), NA)
    
    out <- rbind(out, data.frame(
      Score = s, Beta_unadj = b0, P_unadj = p0,
      Beta_adj = b1, P_adj = p1,
      Attenuation = atten, N = nrow(dat_cc), Note = ""
    ))
  }
  out$P_unadj_FDR <- p.adjust(out$P_unadj, method = "BH")
  out$P_adj_FDR   <- p.adjust(out$P_adj,   method = "BH")
  out <- out[order(out$P_adj, out$P_unadj, na.last = TRUE), ]
  out
}

res_ssg <- run_gee(dat_ssg_cc, score_cols_ssg, "ssGSEA")
res_est <- run_gee(dat_est_cc, score_cols_est, "ESTIMATE")

write.csv(res_ssg, "GEE_ssgsea_adjusted.csv", row.names = FALSE)
write.csv(res_est, "GEE_estimate_adjusted.csv", row.names = FALSE)

cat("\n===== ssGSEA GEE 完成：GEE_ssgsea_adjusted.csv =====\n")
cat("调整后 FDR<0.05 的条目数：", sum(res_ssg$P_adj_FDR < 0.05, na.rm = TRUE),
    " / ", nrow(res_ssg), "\n")
#调整后 FDR<0.05 的条目数： 27  /  29 

cat("\n===== ESTIMATE GEE 完成：GEE_estimate_adjusted.csv =====\n")
print(res_est)
'''
          Score Beta_unadj      P_unadj Beta_adj        P_adj Attenuation   N Note  P_unadj_FDR    P_adj_FDR
3 ESTIMATEScore   549.1821 5.332401e-13 561.6603 9.481305e-14 -0.02272133 860      7.998602e-13 2.844391e-13
1  StromalScore   264.7334 2.731149e-14 261.5154 2.312595e-13  0.01215591 860      8.193446e-14 3.468892e-13
2   ImmuneScore   284.3264 3.210022e-10 302.5985 1.550060e-11 -0.06426449 860      3.210022e-10 1.550060e-11
'''




rm(list=ls())

library(ggplot2)
library(ggpubr)

#CIBERSORT
#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\03immue\\04cibersort'
setwd(path)

source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

#plot
rt = read.table("CIBERSORT-Results.txt",sep="\t",header=T,check.names=F, row.names = 1)

#data=rt[rt[,"P-value"]<0.05, ]
colnames(rt)
#rt<-rt[,-c(1,2,5,6,8,9,12,13,15,16,17,19,20,21,22)]
data<-rt
data = data[ , 1:(ncol(data)-3)]


Type = read.table("cluster.txt",sep="\t",check.names=F,header=T)
rownames(Type)<-Type$id
Type<-Type[,-1]
#names<-rownames(Type)
Type[which(Type$group == 1),]$group<-'Subtype1'
Type[which(Type$group == 2),]$group<-'Subtype2'

colnames(Type)<-c("Subtype","cluster")

Type$id = rownames(Type)
data$id = rownames(data)

data_1 = merge.data.frame(data,Type, by= "id")

data_1 = data_1[order(data_1$Subtype, decreasing = T), ]

#------------------------------
#all
outTab = data.frame()
for(i in colnames(data_1[ , 2:(ncol(data_1)-2)])){
  rt1 = data_1[ , c(i, "Subtype")]
  colnames(rt1) = c("expression","Subtype")
  test = wilcox.test(expression ~ Subtype, data = rt1,exact=FALSE)
  pValue = test$p.value
  outTab = rbind(outTab, cbind(rt1, gene=i))
}

write.table(outTab,file="all_data.txt",sep="\t",row.names=F,quote=F)

#plot all
data=read.table("all_data.txt",sep="\t",header=T,check.names=F)
data[which(data$Subtype == 'Subtype1'),'Subtype']<-'Group1'
data[which(data$Subtype == 'Subtype2'),'Subtype']<-'Group2'
colnames(data)[2]<-'Group'

data$Group=factor(data$Group, levels=c("Group1", "Group2"))
data$expression<-log2(data$expression)

p=ggboxplot(data, x="gene", y="expression", fill = "Group",
            ylab="Fraction",
            xlab="",
            palette = c('lightpink','cadetblue1'))
p=p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 5)+
  theme(text = element_text(size = 12))
p

ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 7, height = 6)


###############################################################################
## 项目：UC 亚型免疫差异 GEE 分析
## 输入：GSE193677.csv、CIBERSORT-Results.txt、NMF_result.csv
## 输出：clinical_metadata.txt、GEE_immune_adjusted.csv
###############################################################################

rm(list = ls())
setwd('D:\\LZ-25234\\UC\\12Subtype\\03immue\\04cibersort')

## 依赖包
if (!requireNamespace("geepack", quietly = TRUE)) install.packages("geepack")
library(geepack)

options(stringsAsFactors = FALSE)

###############################################################################
## Part 1. 从 GSE193677.csv 生成 clinical_metadata.txt
###############################################################################

geo_file <- "GSE193677.csv"   # 若不在当前目录，改成绝对路径

raw <- read.csv(geo_file,
                check.names = FALSE,
                stringsAsFactors = FALSE,
                na.strings = c("NA", ""))

## 关键列检查
need_cols <- c(
  "geo_accession", "title",
  "ibd_disease:ch1", "diseasetypere:ch1", "regionre:ch1", "typere:ch1",
  "ibd_clinicianmeasure_inactive_active:ch1",
  "ibd_endoseverity_4levels:ch1",
  "ibdmesuc_mayo_score:ch1", "ibdsescd_totalsescd:ch1",
  "harveybradshawindex_hbi_score:ch1", "colitisactivityindex_sccai:ch1",
  "nancyindex:ch1", "ghas_sum7:ch1",
  "endoremiss:ch1", "historemiss:ch1",
  "demographics_gender:ch1", "study_eligibility_age_at_endo:ch1",
  "crp_jjmgl_log2:ch1", "log2_fecalcalpro_mgperg:ch1"
)
missing_cols <- setdiff(need_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("GSE193677.csv 缺少列：", paste(missing_cols, collapse = ", "))
}

## 整理临床表
clin_all <- data.frame(
  id                  = raw$geo_accession,
  Patient_ID          = sub("_Biopsy_.*$", "", raw$title),
  Disease             = raw[["ibd_disease:ch1"]],
  DiseaseType         = raw[["diseasetypere:ch1"]],
  BiopsySite          = raw[["regionre:ch1"]],
  Activity            = raw[["ibd_clinicianmeasure_inactive_active:ch1"]],
  Endoscopy4          = raw[["ibd_endoseverity_4levels:ch1"]],
  Mayo                = suppressWarnings(as.numeric(raw[["ibdmesuc_mayo_score:ch1"]])),
  SESCD               = suppressWarnings(as.numeric(raw[["ibdsescd_totalsescd:ch1"]])),
  HBI                 = suppressWarnings(as.numeric(raw[["harveybradshawindex_hbi_score:ch1"]])),
  SCCAI               = suppressWarnings(as.numeric(raw[["colitisactivityindex_sccai:ch1"]])),
  Nancy               = suppressWarnings(as.numeric(raw[["nancyindex:ch1"]])),
  GHAS                = suppressWarnings(as.numeric(raw[["ghas_sum7:ch1"]])),
  Inflammation        = raw[["typere:ch1"]],
  EndoscopicRemission = raw[["endoremiss:ch1"]],
  HistologicRemission = raw[["historemiss:ch1"]],
  Gender              = raw[["demographics_gender:ch1"]],
  Age                 = suppressWarnings(as.numeric(raw[["study_eligibility_age_at_endo:ch1"]])),
  CRP_log2            = suppressWarnings(as.numeric(raw[["crp_jjmgl_log2:ch1"]])),
  FecalCal_log2       = suppressWarnings(as.numeric(raw[["log2_fecalcalpro_mgperg:ch1"]])),
  stringsAsFactors    = FALSE
)

## 因子化：Inflammation 用 typere 的 I/NonI
clin_all$Inflammation <- factor(clin_all$Inflammation,
                                levels = c("NonI", "I"),
                                labels = c("NonInflamed", "Inflamed"))
clin_all$Activity   <- factor(clin_all$Activity, levels = c("Inactive", "Active"))
clin_all$BiopsySite <- factor(clin_all$BiopsySite)
clin_all$Patient_ID <- factor(clin_all$Patient_ID)
clin_all$Disease    <- factor(clin_all$Disease)

## 仅保留 UC（若你的 NMF 分组包含 CD/Control，请改成 FALSE 并自行分疾病建模）
uc_only <- TRUE
clin <- if (uc_only) subset(clin_all, Disease == "UC") else clin_all
clin$Patient_ID <- droplevels(clin$Patient_ID)

## Endoscopy：UC 主用 Mayo 内镜评分 0-3
clin$Endoscopy <- clin$Mayo
if (mean(is.na(clin$Endoscopy)) > 0.5) {
  message("提示：Mayo 缺失超过 50%，建议改用 Endoscopy4 或检查原始临床表。")
}
clin$Endoscopy <- factor(clin$Endoscopy,
                         levels = sort(unique(clin$Endoscopy[!is.na(clin$Endoscopy)])))

## 输出文件（Tab 分隔，含 id 列，无行名）
clin_out <- clin[, c(
  "id", "Patient_ID", "Activity", "Endoscopy", "Inflammation", "BiopsySite",
  "Disease", "Gender", "Age",
  "Mayo", "Nancy", "GHAS",
  "EndoscopicRemission", "HistologicRemission",
  "CRP_log2", "FecalCal_log2"
)]

if (anyDuplicated(clin_out$id)) stop("clinical_metadata 中 id 有重复")

write.table(clin_out, "clinical_metadata.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## 临床变量分布与缺失检查
cat("\n===== clinical_metadata.txt 生成完成 =====\n")
cat("样本数：", nrow(clin_out), "，患者数：", length(unique(clin_out$Patient_ID)), "\n")
#样本数： 872 ，患者数： 421 

cat("\nGroup 相关变量缺失情况：\n")
print(colSums(is.na(clin_out[, c("Activity", "Endoscopy", "Inflammation", "BiopsySite", "Patient_ID")])))

cat("\nActivity 分布：\n");       print(table(clin_out$Activity, useNA = "ifany"))
'''
Activity 分布：

Inactive   Active 
     746      126 
'''
cat("\nEndoscopy 分布：\n");      print(table(clin_out$Endoscopy, useNA = "ifany"))
'''
Endoscopy 分布：

   0    1    2    3 <NA> 
 337  315  143   72    5 
'''
cat("\nInflammation 分布：\n");   print(table(clin_out$Inflammation, useNA = "ifany"))
'''
Inflammation 分布：

NonInflamed    Inflamed 
        579         293 
'''
cat("\nBiopsySite 分布：\n");     print(table(clin_out$BiopsySite, useNA = "ifany"))
'''
BiopsySite 分布：

     Cecum      Ileum  LeftColon     Rectum RightColon    Sigmoid Transverse 
        66        188         90        311         95         89         33 
'''

###############################################################################
## Part 2. 读入 CIBERSORT、NMF 分组，并与临床表合并
###############################################################################

## 2.1 CIBERSORT
cib <- read.table("CIBERSORT-Results.txt",
                  sep = "\t", header = TRUE, row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)

## 按列名删掉 CIBERSORT 统计列，避免硬删位置
drop_cols <- intersect(c("P-value", "Correlation", "RMSE"), colnames(cib))
if (length(drop_cols) > 0) cib <- cib[, setdiff(colnames(cib), drop_cols), drop = FALSE]

cib$id <- rownames(cib)
cell_types <- setdiff(colnames(cib), "id")
cat("\nCIBERSORT 细胞类型数：", length(cell_types), "\n")
#CIBERSORT 细胞类型数： 22 
cat("CIBERSORT id 示例：", paste(head(cib$id, 3), collapse = ", "), "\n")

## 检查全 0 细胞类型
zero_cells <- cell_types[colSums(cib[, cell_types, drop = FALSE]) == 0]
if (length(zero_cells) > 0) {
  warning("以下细胞类型全为 0，将从分析中剔除：", paste(zero_cells, collapse = ", "))
  cell_types <- setdiff(cell_types, zero_cells)
}

## 2.2 NMF 分组
grp <- read.csv("NMF_result.csv",
                check.names = FALSE, stringsAsFactors = FALSE)
colnames(grp) <- c('id', 'Group')
grp[which(grp$Group == 1), 'Group'] <- 'Group1'
grp[which(grp$Group == 2), 'Group'] <- 'Group2'

if (!all(c("id", "Group") %in% colnames(grp))) {
  stop("NMF_result.csv 必须包含 id 和 Group 两列")
}
grp$id    <- as.character(grp$id)
grp$Group <- factor(grp$Group, levels = c("Group1", "Group2"))
if (any(is.na(grp$Group))) stop("NMF_result.csv 的 Group 只能包含 Group1 / Group2")
cat("NMF 分组样本数：", nrow(grp), "\n")
#NMF 分组样本数： 865
cat("NMF id 示例：", paste(head(grp$id, 3), collapse = ", "), "\n")

## 2.3 读临床表（保留 id 列，不要 row.names=1）
clin <- read.table("clinical_metadata.txt",
                   sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
clin$id <- as.character(clin$id)

## 修改：read.table 读回后分类变量为 character，统一还原为 factor
clin$Activity     <- factor(clin$Activity, levels = c("Inactive", "Active"))
clin$Inflammation <- factor(clin$Inflammation, levels = c("NonInflamed", "Inflamed"))
clin$BiopsySite   <- factor(clin$BiopsySite)
clin$Patient_ID   <- factor(clin$Patient_ID)
clin$Endoscopy    <- factor(clin$Endoscopy,
                            levels = sort(unique(clin$Endoscopy[!is.na(clin$Endoscopy)])))

## 2.4 ID 对齐检查
cat("\n===== ID 对齐检查 =====\n")
cat("CIBERSORT ∩ NMF：", length(intersect(cib$id, grp$id)), "\n")
#CIBERSORT ∩ NMF： 865 
cat("CIBERSORT ∩ clinical：", length(intersect(cib$id, clin$id)), "\n")
#CIBERSORT ∩ clinical： 865
cat("NMF ∩ clinical：", length(intersect(grp$id, clin$id)), "\n")
#NMF ∩ clinical： 865 

if (length(intersect(cib$id, grp$id)) == 0) {
  stop("CIBERSORT 与 NMF_result.csv 的 id 无交集，请检查 id 格式（GSM 号 vs 样本名）")
}
if (length(intersect(cib$id, clin$id)) == 0) {
  stop("CIBERSORT 与 clinical_metadata.txt 的 id 无交集，请检查 id 格式")
}

## 2.5 合并
dat <- merge(cib, grp,  by = "id")
dat <- merge(dat, clin, by = "id")
cat("\n合并后样本数：", nrow(dat), "\n")
#合并后样本数： 865 

## 2.6 完整案例集（调整模型使用）
gee_vars <- c("Group", "Activity", "Endoscopy", "Inflammation", "BiopsySite", "Patient_ID")
dat_cc <- dat[complete.cases(dat[, gee_vars, drop = FALSE]), ]

## 修改：每列先 factor 再 droplevels，避免 character 报错
dat_cc$Group        <- droplevels(factor(dat_cc$Group, levels = c("Group1", "Group2")))
dat_cc$Patient_ID   <- droplevels(factor(dat_cc$Patient_ID))
dat_cc$Activity     <- droplevels(factor(dat_cc$Activity, levels = c("Inactive", "Active")))
dat_cc$Endoscopy    <- droplevels(factor(dat_cc$Endoscopy,
                                         levels = sort(unique(dat_cc$Endoscopy[!is.na(dat_cc$Endoscopy)]))))
dat_cc$Inflammation <- droplevels(factor(dat_cc$Inflammation, levels = c("NonInflamed", "Inflamed")))
dat_cc$BiopsySite   <- droplevels(factor(dat_cc$BiopsySite))

cat("GEE 完整案例样本数：", nrow(dat_cc),
    "，患者数：", length(unique(dat_cc$Patient_ID)), "\n")
#GEE 完整案例样本数： 860 ，患者数： 416

cat("完整案例中 Group 分布：\n"); print(table(dat_cc$Group, useNA = "ifany"))
'''
完整案例中 Group 分布：

Group1 Group2 
   472    388 
'''
if (nrow(dat_cc) == 0) stop("完整案例集为 0，请检查 Activity/Endoscopy/Inflammation/BiopsySite 缺失")

## 修改：Group 必须两个水平都在，否则 GEE 没有 GroupGroup2 系数
if (nlevels(dat_cc$Group) < 2) {
  stop("完整案例集中 Group 只有一个水平，无法比较 Group2 vs Group1")
}

## 修改：分类协变量至少两个水平，否则该协变量无法进入模型
for (v in c("Activity", "Endoscopy", "Inflammation", "BiopsySite")) {
  if (nlevels(dat_cc[[v]]) < 2) {
    warning("完整案例集中 ", v, " 只有一个水平，调整模型可能失败")
  }
}

###############################################################################
## Part 3. GEE：未调整 vs 调整模型
###############################################################################

gee_res <- data.frame()

for (cell in cell_types) {
  f0 <- as.formula(paste0("`", cell, "` ~ Group"))
  f1 <- as.formula(paste0("`", cell, "` ~ Group + Activity + Endoscopy + Inflammation + BiopsySite"))
  
  m0 <- tryCatch(
    geeglm(f0, data = dat_cc, id = Patient_ID, family = gaussian, corstr = "exchangeable"),
    error = function(e) NULL
  )
  m1 <- tryCatch(
    geeglm(f1, data = dat_cc, id = Patient_ID, family = gaussian, corstr = "exchangeable"),
    error = function(e) NULL
  )
  
  if (is.null(m0) || is.null(m1)) {
    gee_res <- rbind(gee_res, data.frame(
      CellType = cell, Beta_unadj = NA, P_unadj = NA,
      Beta_adj = NA, P_adj = NA, Attenuation = NA, N = nrow(dat_cc),
      Note = "模型拟合失败"
    ))
    next
  }
  
  coef0 <- summary(m0)$coefficients
  coef1 <- summary(m1)$coefficients
  
  ## 防止 GroupGroup2 不存在（例如某个水平缺失）
  if (!"GroupGroup2" %in% rownames(coef0) || !"GroupGroup2" %in% rownames(coef1)) {
    gee_res <- rbind(gee_res, data.frame(
      CellType = cell, Beta_unadj = NA, P_unadj = NA,
      Beta_adj = NA, P_adj = NA, Attenuation = NA, N = nrow(dat_cc),
      Note = "GroupGroup2 系数缺失"
    ))
    next
  }
  
  b0 <- coef0["GroupGroup2", "Estimate"]
  p0 <- coef0["GroupGroup2", "Pr(>|W|)"]
  b1 <- coef1["GroupGroup2", "Estimate"]
  p1 <- coef1["GroupGroup2", "Pr(>|W|)"]
  
  atten <- ifelse(abs(b0) > 1e-4, 1 - abs(b1) / abs(b0), NA)
  
  gee_res <- rbind(gee_res, data.frame(
    CellType = cell,
    Beta_unadj = b0, P_unadj = p0,
    Beta_adj   = b1, P_adj   = p1,
    Attenuation = atten,
    N = nrow(dat_cc),
    Note = ""
  ))
}

## 整理输出
gee_res$P_unadj_FDR <- p.adjust(gee_res$P_unadj, method = "BH")
gee_res$P_adj_FDR   <- p.adjust(gee_res$P_adj,   method = "BH")

## 按调整后 P 值排序
gee_res <- gee_res[order(gee_res$P_adj, gee_res$P_unadj, na.last = TRUE), ]

write.csv(gee_res, "GEE_immune_adjusted.csv", row.names = FALSE)

cat("\n===== GEE 分析完成 =====\n")
cat("输出文件：GEE_immune_adjusted.csv\n")
cat("\nMacrophage / Mast / Dendritic / Plasma 相关结果：\n")
print(gee_res[grepl("Macrophage|Mast|Dendritic|Plasma", gee_res$CellType), ])

###############################################################################
## Part 4. 可选：敏感性分析（Endoscopy 用四分类，Inflammation 用 Nancy）
###############################################################################

## 如果 Mayo 缺失较多，可把 Endoscopy 换成 Endoscopy4 再跑一次
## 如果 Nancy 缺失可接受，可把 Inflammation 换成 Nancy 二分类（如 Nancy>=2 vs <2）

# dat_sens <- dat
# dat_sens$Endoscopy <- factor(dat_sens$Endoscopy4)
# dat_sens$Inflammation <- factor(ifelse(dat_sens$Nancy >= 2, "Active", "Quiescent"))
# dat_sens$Patient_ID <- factor(dat_sens$Patient_ID)
# dat_sens <- dat_sens[complete.cases(dat_sens[, gee_vars, drop = FALSE]), ]
# 然后复用 Part 3 循环，输出 GEE_immune_sensitivity.csv





rm(list=ls())

library(limma)
library(edgeR)
library(doBy)


path<-'D:\\LZ-25234\\UC\\12Subtype\\04diff'
setwd(path)

data<-read.table('normalize.txt',header = T,row.names=1,sep='\t')
Exp<-data

group_list<-read.csv('NMF_result.csv',header=T,sep=',')
table(group_list$group)
# 1  2 
#474 391 

group_list[which(group_list$group == 1),]$group<-'group1'
group_list[which(group_list$group == 2),]$group<-'group2'
group_list<-group_list[order(group_list$group),]
head(group_list)

cols<-group_list$X
Exp<-Exp[,cols]
rownames(Exp)<-toupper(rownames(Exp))

group_list<-group_list$group
group<-as.factor(group_list)
design<-model.matrix(~0+group)
rownames(design)<-colnames(Exp)
head(design)

#分组文件
ALL_group<-data.frame('Sample'=colnames(Exp),'Type'=group_list)
dim(ALL_group)
#[1]  865 2

write.table(ALL_group,'group_list.txt',row.names=F,col.names=T,quote=F,sep='\t')


dim(Exp)
table(group)
#group1 group2 
#   474    391

#无需进行voom转换
#线性模型拟合
exprSet<-Exp
fit<-lmFit(exprSet, design)

contrast.matrix <- makeContrasts(groupgroup2 - groupgroup1, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, adjust.method="BH", number=Inf)
dim(deg)  #[1] 39374     6

#保存所有差异分析的结果
write.table(deg,'all.xls',col.names = T,row.names = T,sep = '\t')

#Filtering thresholds
adjP=0.05
fc=1

output<-deg
#result1-fc=0.584
outdiff<-output[(output$adj.P.Val < adjP & (output$logFC>fc | output$logFC<(-fc))),]
#outdiff<-output[(output$P.Value < adjP),]
dim(outdiff)
#[1] 1083   6

#基因上下调
a<-outdiff
a$Group<-'stable'
a[which(a$logFC > 0),'stable']<-'up'
a[which(a$logFC < 0),'stable']<-'down'
table(a$stable)
#down   up 
# 460   623 

head(outdiff,10)

write.table(outdiff,'outdiff_P0.05.txt',row.names = T,sep = '\t',quote = F)

#差异基因表达谱
Normalize_counts<-read.table("normalize.txt", sep = "\t", header = T)
outdiff<-read.table('outdiff_P0.05.txt',header = T,sep = '\t')
data<-Normalize_counts
heatmap2=data[rownames(outdiff),]
dim(heatmap2) #[1] 1083    1320

#save
write.table(heatmap2,file="diffGeneExp_p0.05.txt",sep="\t",col.names=T,row.names=T,quote=F)

###火山图绘制
res<-deg
logFC<-as.numeric(res$logFC)
FDR<-res$adj.P.Val


#p<0.05 logfc>1
plot(logFC, -1*log10(FDR),xlim=c(-2,2),
     col=ifelse(FDR<=0.05 & abs(logFC)>1, "red", "black"),
     xlab="logFC", ylab="-1*log10(adj.P.Val)", main="Volcano plot", pch=".")


pdf(file="Volcano.pdf", width = 4, height = 4)
plot(logFC, -1*log10(FDR), 
     xlim=c(-3,3),
     col=ifelse(FDR <= 0.05 & abs(logFC) > 1,
                ifelse(logFC > 1, "red", "blue"),
                "black"),
     xlab="logFC", ylab="-1*log10(adj.P.Val)", main="Volcano plot", pch=".")
dev.off()




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\05clusterprofiler'
setwd(path)

#BiocManager::install('GOplot')
library(GOplot)#弦图，弦表图，系统聚类图
library(clusterProfiler)
library(stringr)

#data
#info<-read.table('diff_gse193677_logfc_0.585.xls',header=T,row.names=1,sep='\t')
info<-read.table('outdiff_P0.05.txt',header=T,sep='\t')
genes<-rownames(info)
length(genes)
#[1] 1083

#指定富集分析的物种库
GO_database<-'org.Hs.eg.db' #GO分析指定物种
KEGG_database<-'hsa' #KEGG分析指定物种

#gene ID转换
gene<-bitr(genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
dim(gene)
#[1] 1048    2


##GO分析 这里只做GO
GO<-enrichGO(gene$ENTREZID,#GO富集分析
  OrgDb = GO_database,
  keyType = "ENTREZID",#设定读取的gene ID类型
  ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
  pvalueCutoff = 0.05,#设定p值阈值
  qvalueCutoff = 0.05,#设定q值阈值
  readable = T)

#查看一下GO的情况
result<-GO@result
dim(result) #[1] 458   13
table(result$ONTOLOGY)
# BP  CC  MF 
#306  40  112 

#导出结果
write.csv(result,'go_result.csv',quote = F)

#这里确认了BP从1:306,CC从307:346,MF从347:458

##这里可能要对数据做一些处理方便定位
genedata<-data.frame()
genedata<-data.frame('ID'=genes,'logFC'=info$logFC)

##弦图绘制
GOplotIn_BP<-GO[1:5,c(2,3,10,12)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[307:311,c(2,3,10,12)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[347:351,c(2,3,10,12)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列

GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')

names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')

GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"

circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 

chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 

'''
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)
'''

##聚类图
chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
'''
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 
'''

#保存图片
tiff(file="GO_BP.tiff",width = 15,height = 10,units = "in",bg = "white",res = 500)
chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
dev.off()

pdf(file = "GO_BP.pdf",height = 10,width = 15)
chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
dev.off()

'''
tiff(file="GO_CC.tiff",width = 15,height = 10,units = "in",bg = "white",res = 500)
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
dev.off()

pdf(file = "GO_CC.pdf",height = 10,width = 15)
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
dev.off()

tiff(file="GO_MF.tiff",width = 15,height = 10,units = "in",bg = "white",res = 500)
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 
dev.off()

pdf(file = "GO_MF.pdf",height = 10,width = 15)
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 
dev.off()
'''

#圈图
GOCircle(circ_BP)
#circ_BP<-circ_BP[!duplicated(circ_BP),]

#保存图片

pdf(file = "GO_BP.pdf",height = 6,width = 14)
GOCircle(circ_BP)
dev.off()





rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\12Subtype\\05clusterprofiler'
setwd(path)

#data
diff<-read.table('outdiff_P0.05.txt',header=T,sep='\t')
dim(diff)
#[1] 1083   6

#这里没有区分上下调组
high<-diff

#通路分析
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)

library(GOplot)#弦图，弦表图，系统聚类图
library(stringr)


#查看支持转换的类型
keytypes(org.Hs.eg.db)

#指定富集分析的物种库
GO_database<-'org.Hs.eg.db' #GO分析指定物种
KEGG_database<-'hsa' #KEGG分析指定物种

#gene ID转换
genes<-rownames(high)
gene<-bitr(genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
dim(gene)
#[1] 1048    2


##GSEA的KEGG分析
info<-diff
info<-info[,c(1,4,5)]
info$Sample<-rownames(info)
colnames(info)[1:3]<-c('logFC','pValue','fdr')
info<-info[,c('Sample','logFC','pValue','fdr')]
names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange

GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)#GSEA富集分析


#这里保存结果
write.csv(GSEA_KEGG@result,'KEGG_result.csv',quote = F)

#plot  这里KEGG富集到的结果只有1个
ridgeplot(GSEA_KEGG,10)

'''
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,1:4)#30是根据ridgeplot中有30个富集通路得到的
'''
pdf(file = "KEGG.pdf",height = 4,width = 6)
ridgeplot(GSEA_KEGG,10)
dev.off()

#下面绘图 调整坐标轴大小
library(enrichplot)
library(ggplot2)

p <- ridgeplot(GSEA_KEGG, 10) + 
  theme(
    axis.text = element_text(size = 8),      # 坐标轴刻度文字大小
    axis.title = element_text(size = 16),     # 坐标轴标题大小
    plot.title = element_text(size = 15),     # 图标题大小
    legend.text = element_text(size = 12),    # 图例文字大小
    legend.title = element_text(size = 14)    # 图例标题大小
  )

print(p)

pdf(file = "KEGG.pdf",height = 6,width = 6)
p
dev.off()




rm(list = ls())

# ============================================================
# UC subtype clinical associations and percentage bar plots
# Primary GEE working correlation: independence
# Patient-level clustering with robust sandwich standard errors
# ============================================================

# 0. Packages
packages <- c("dplyr", "tidyr", "ggplot2", "geepack", "scales")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(geepack)
library(scales)

# 1. Working directory and input files
# Set this to the folder containing GSE193677.csv and NMF_result.csv
path <- "D:/LZ-25234/UC/12Subtype/06barplot"
setwd(path)

metadata <- read.csv(
  "GSE193677.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cluster <- read.csv(
  "NMF_result.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 2. Merge metadata and subtype assignments
colnames(cluster)[1] <- "geo_accession"

metadata$geo_accession <- as.character(metadata$geo_accession)
cluster$geo_accession <- as.character(cluster$geo_accession)

dat <- metadata %>%
  inner_join(cluster, by = "geo_accession")

names(dat)[1] <- "sample_index"

stopifnot(nrow(dat) == 865)
stopifnot(sum(duplicated(dat$geo_accession)) == 0)

# 3. Extract patient ID and clinical variables
dat <- dat %>%
  mutate(
    title_prefix = sub(",.*$", "", title),
    patient_id = sub("_Biopsy_[0-9]+$", "", title_prefix),

    group = factor(
      as.character(group),
      levels = c("1", "2"),
      labels = c("Group 1", "Group 2")
    ),

    activity = trimws(as.character(
      `ibd_clinicianmeasure_inactive_active:ch1`
    )),
    endo_severity = trimws(as.character(
      `ibd_endoseverity_4levels:ch1`
    )),
    tissue_status = trimws(as.character(`typere:ch1`)),
    region = trimws(as.character(`regionre:ch1`)),

    group_binary = ifelse(group == "Group 2", 1, 0),

    active_binary = case_when(
      activity == "Active" ~ 1,
      activity == "Inactive" ~ 0,
      TRUE ~ NA_real_
    ),

    endo_binary = case_when(
      endo_severity %in% c("Moderate", "Severe") ~ 1,
      endo_severity %in% c("Inactive", "Mild") ~ 0,
      TRUE ~ NA_real_
    ),

    inflamed_binary = case_when(
      tissue_status == "I" ~ 1,
      tissue_status == "NonI" ~ 0,
      TRUE ~ NA_real_
    )
  )

dat$group <- factor(dat$group, levels = c("Group 1", "Group 2"))
dat$region <- factor(dat$region)

cat("Merged samples:", nrow(dat), "\n")
#Merged samples: 865 
cat("Unique patients:", dplyr::n_distinct(dat$patient_id), "\n")
#Unique patients: 419 
print(table(dat$group))
#Group 1 Group 2 
# 474     391 

# 4. Clinical GEE function
# Independence is the working correlation structure.
# Patient clustering is retained through id_cluster.
run_gee_ind <- function(data, outcome, adjust_region = FALSE) {

  if (adjust_region) {
    formula <- as.formula(paste0(outcome, " ~ group + region"))
  } else {
    formula <- as.formula(paste0(outcome, " ~ group"))
  }

  required_vars <- c(outcome, "group", "patient_id")
  if (adjust_region) required_vars <- c(required_vars, "region")

  tmp <- data %>%
    select(all_of(required_vars)) %>%
    filter(complete.cases(.)) %>%
    mutate(id_cluster = as.integer(factor(patient_id))) %>%
    arrange(id_cluster)

  fit <- geeglm(
    formula,
    id = id_cluster,
    data = tmp,
    family = binomial(link = "logit"),
    corstr = "independence",
    std.err = "san.se"
  )

  beta <- coef(fit)
  V <- vcov(fit)
  term <- "groupGroup 2"

  estimate <- beta[term]
  variance <- V[term, term]

  if (!is.finite(variance) || variance <= 0) {
    stop("Invalid robust variance for outcome: ", outcome)
  }

  se <- sqrt(variance)

  result <- data.frame(
    Outcome = outcome,
    Model = ifelse(adjust_region, "Adjusted for region", "Unadjusted"),
    Working_correlation = "independence",
    N = nrow(tmp),
    Patients = dplyr::n_distinct(tmp$id_cluster),
    OR = exp(estimate),
    Lower_95CI = exp(estimate - 1.96 * se),
    Upper_95CI = exp(estimate + 1.96 * se),
    P_value = 2 * pnorm(-abs(estimate / se))
  )

  list(fit = fit, result = result)
}

# 5. Clinical association analyses
ind_activity_1 <- run_gee_ind(dat, "active_binary", FALSE)
ind_activity_2 <- run_gee_ind(dat, "active_binary", TRUE)

ind_endo_1 <- run_gee_ind(dat, "endo_binary", FALSE)
ind_endo_2 <- run_gee_ind(dat, "endo_binary", TRUE)

ind_inflamed_1 <- run_gee_ind(dat, "inflamed_binary", FALSE)
ind_inflamed_2 <- run_gee_ind(dat, "inflamed_binary", TRUE)

clinical_results <- bind_rows(
  ind_activity_1$result,
  ind_activity_2$result,
  ind_endo_1$result,
  ind_endo_2$result,
  ind_inflamed_1$result,
  ind_inflamed_2$result
) %>%
  group_by(Model) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup()

print(clinical_results)

write.csv(
  clinical_results,
  "Clinical_GEE_independence_results.csv",
  row.names = FALSE
)

# 6. Subtype versus biopsy region: independence GEE
region_data <- dat %>%
  filter(
    !is.na(group_binary),
    !is.na(region),
    !is.na(patient_id)
  ) %>%
  mutate(id_cluster = as.integer(factor(patient_id))) %>%
  arrange(id_cluster)

region_fit <- geeglm(
  group_binary ~ region,
  id = id_cluster,
  data = region_data,
  family = binomial(link = "logit"),
  corstr = "independence",
  std.err = "san.se"
)

print(summary(region_fit))

# Overall robust Wald test for all region coefficients
beta <- coef(region_fit)
V <- vcov(region_fit)

region_terms <- grep("^region", names(beta), value = TRUE)
beta_region <- beta[region_terms]
V_region <- V[region_terms, region_terms, drop = FALSE]

if (any(!is.finite(V_region)) ||
    any(!is.finite(beta_region)) ||
    qr(V_region)$rank < length(region_terms)) {
  stop("Region robust covariance matrix is singular or invalid.")
}

wald_stat <- as.numeric(
  t(beta_region) %*% solve(V_region, beta_region)
)

df_region <- length(region_terms)
p_region <- pchisq(wald_stat, df = df_region, lower.tail = FALSE)

region_result <- data.frame(
  N = nrow(region_data),
  Patients = dplyr::n_distinct(region_data$id_cluster),
  Working_correlation = "independence",
  Wald_statistic = wald_stat,
  df = df_region,
  P_value = p_region
)

print(region_result)
'''
    N Patients Working_correlation Wald_statistic df   P_value
1 865      419        independence          73.42  6 8.123e-14
'''

write.csv(
  region_result,
  "Subtype_region_GEE_independence_results.csv",
  row.names = FALSE
)

# 7. Percentage stacked bar plots
# Percentages are calculated within each subtype group.
# Labels are placed in the middle of each stacked segment.
plot_clinical <- function(data, variable, plot_title, filename) {

  plot_data <- data %>%
    filter(
      !is.na(.data[[variable]]),
      !is.na(group)
    ) %>%
    count(group, category = .data[[variable]], name = "n") %>%
    group_by(group) %>%
    mutate(
      percentage = n / sum(n),
      percentage_label = scales::percent(percentage, accuracy = 0.1)
    ) %>%
    ungroup()

  p <- ggplot(
    plot_data,
    aes(x = group, y = percentage, fill = category)
  ) +
    geom_col(
      width = 0.65,
      color = "black",
      linewidth = 0.3
    ) +
    geom_text(
      aes(label = percentage_label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      color = "black"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Percentage",
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "right"
    )

  ggsave(
    filename,
    p,
    width = 5.5,
    height = 5,
    dpi = 300
  )

  # Also save a high-resolution PNG version
  png_filename <- sub("\\.pdf$", ".png", filename)
  ggsave(
    png_filename,
    p,
    width = 5.5,
    height = 5,
    dpi = 600
  )

  return(p)
}

# Disease activity
p1 <- plot_clinical(
  dat,
  "activity",
  "Disease activity",
  "Subtype_disease_activity_percentage.pdf"
)

# Endoscopic severity
p2 <- plot_clinical(
  dat,
  "endo_severity",
  "Endoscopic severity",
  "Subtype_endoscopic_severity_percentage.pdf"
)

# Tissue inflammation
p3 <- plot_clinical(
  dat,
  "tissue_status",
  "Tissue inflammation",
  "Subtype_tissue_inflammation_percentage.pdf"
)

# Biopsy location
p4 <- plot_clinical(
  dat,
  "region",
  "Biopsy location",
  "Subtype_biopsy_location_percentage.pdf"
)

print(p1)
print(p2)
print(p3)
print(p4)

# 8. Export descriptive percentages for checking/reporting
make_percentage_table <- function(data, variable) {
  data %>%
    filter(!is.na(.data[[variable]]), !is.na(group)) %>%
    count(group, category = .data[[variable]], name = "n") %>%
    group_by(group) %>%
    mutate(
      denominator = sum(n),
      percentage = 100 * n / denominator
    ) %>%
    ungroup()
}

percentage_tables <- bind_rows(
  make_percentage_table(dat, "activity") %>% mutate(Variable = "activity"),
  make_percentage_table(dat, "endo_severity") %>% mutate(Variable = "endo_severity"),
  make_percentage_table(dat, "tissue_status") %>% mutate(Variable = "tissue_status"),
  make_percentage_table(dat, "region") %>% mutate(Variable = "region")
) %>%
  select(Variable, group, category, n, denominator, percentage)

write.csv(
  percentage_tables,
  "Subtype_clinical_and_region_percentages.csv",
  row.names = FALSE
)

cat("\nAnalysis complete. Output files:\n")
cat("- Clinical_GEE_independence_results.csv\n")
cat("- Subtype_region_GEE_independence_results.csv\n")
cat("- Subtype_clinical_and_region_percentages.csv\n")
cat("- Four percentage-labelled PDF and PNG plots\n")

