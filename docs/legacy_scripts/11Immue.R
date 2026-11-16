

rm(list=ls())

#TCGA免疫微环境分析和免疫治疗反应性预测

library(GSVA)
library(limma)
library(GSEABase)


#path
path<-'D:\\LZ-25234\\UC\\11Immue\\01data'
setwd(path)

#data
rt<-read.table('normalize.txt',sep = '\t',header = T,row.names = 1)
Risk_sample<-read.table('group.txt',header = T,sep = '\t',check.names=F)
dim(rt);dim(Risk_sample)
#[1] 39374  1320
#[1] 1320    2

Risk_sampleid<-Risk_sample$Sample
rt<-rt[,Risk_sampleid]
dim(rt) #[1] 39374  1320

##要把基因名换成大写
rt$gene<-rownames(rt)
rt$gene<-toupper(rt$gene)
#处理矩阵
merge_data<-rt
gene<-data.frame(merge_data$gene)
exp<-merge_data[,1:ncol(merge_data)-1]
dim(gene)
dim(exp)
exp[1:4,1:3]

exp<-cbind(gene,exp)
colnames(exp)[1]<-"Symbol"
dim(exp)
#[1] 39374  1321

exp[1:3,1:4]

rt<-exp
rt<-as.matrix(rt)
rownames(rt)<-rt[,1]
exp<-rt[,2:ncol(rt)]
dimnames<-list(rownames(exp),colnames(exp))
exp<-matrix(as.numeric(exp),nrow = nrow(exp), dimnames = dimnames)
dim(exp)  #[1] 39374  1320

exp[1:3,1:3]

#Deduplication
counts<-exp
counts_delrep<-counts
counts_delrep<-avereps(counts_delrep)
Exp<-counts_delrep
dim(Exp)
#[1] 39374  1320

exp<-Exp
dimnames=list(rownames(exp), colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0, ]
dim(mat) #[1] 39374  1320

gmtFile="immune.gmt"

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




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\11Immue\\03estimateHM'
setwd(path)

rt=read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)   #¶ÁÈ¡ÎÄ¼þ
Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=T)

rt=t(rt[,row.names(Type)])
unique(rownames(rt) == rownames(Type))
colnames(rt)

num_immu_cell=c(1,4,6,9,11,12,14,15,17,20,21,22,23,24,25)
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
data$Subtype = factor(data$Subtype,levels=c('Control','CIRI'))
colnames(data)[3]<-'Group'
#data$expression<-log2(data$expression)

p = ggboxplot(data, x="gene", y="expression", fill = "Group",
              ylab = "Fraction",
              xlab = "", palette = c("#4DBBD5FF","#DC0000FF"))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 10)+
  theme(text = element_text(size = 25))
p+coord_flip()

ggsave(filename = "immu_cell_boxplot.pdf", plot = p, dpi = 600, width = 10,
       height = 11)


#plot num_immu_func
data=read.table("data_func.txt",sep="\t",header=T,check.names=F)  
data$Subtype = factor(data$Subtype,levels=c('Control','CIRI'))
colnames(data)[3]<-'Group'

p = ggboxplot(data, x="gene", y="expression", fill = "Group",
              ylab = "Fraction",
              xlab = "", palette = c("#4DBBD5FF","#DC0000FF"))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 10)+
  theme(text = element_text(size = 25))
p

ggsave(filename = "immu_func_boxplot.tiff", plot = p, dpi = 500, width = 10,
       height = 11)
ggsave(filename = "immu_func_boxplot.pdf", plot = p, dpi = 600, width = 10,
       height = 11)



rm(list=ls())

#CIBERSORT
#path
path<-'D:\\LZ-25234\\UC\\11Immue\\04cibersort'
setwd(path)
library(ggpubr)

source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

#plot
rt = read.table("CIBERSORT-Results.txt",sep="\t",header=T,check.names=F, row.names = 1) 
dim(rt) #[1] 1320 25


'''
#data=rt[rt[,"P-value"]<0.05, ]
nonzero_counts <- colSums(rt != 0)
filtered_rt <- rt[, nonzero_counts > 5]

data<-filtered_rt
'''

data<-rt
data = data[ , 1:(ncol(data)-3)]

Type = read.table("cluster.txt",sep="\t",check.names=F,header=T)
rownames(Type)<-Type$Sample
Type<-Type[,-1]
names<-rownames(Type)

colnames(Type)<-c("Subtype","cluster")

Type$id = rownames(Type)
data$id = rownames(data)

data_1 = merge.data.frame(data,Type, by= "id")

data_1 = data_1[order(data_1$Subtype, decreasing = T), ]
colnames(data_1)[24]<-'Group'

#data_1<-data_1[,-8]

#p<0.05
outTab = data.frame()
for(i in colnames(data_1[ , 2:(ncol(data_1)-2)])){
  rt1 = data_1[ , c(i, "Group")]
  colnames(rt1) = c("expression","Group")
  test = wilcox.test(expression ~ Group, data = rt1,exact=FALSE)
  pValue = test$p.value
  if(!is.na(pValue) & pValue < 0.05){
    outTab = rbind(outTab, cbind(rt1, gene=i))
    print(pValue)
  }
}

write.table(outTab,file="P0.05_data.txt",sep="\t",row.names=F,quote=F)


#all
outTab = data.frame()
for(i in colnames(data_1[ , 2:(ncol(data_1)-2)])){
  rt1 = data_1[ , c(i, "Group")]
  colnames(rt1) = c("expression","Group")
  test = wilcox.test(expression ~ Group, data = rt1,exact=FALSE)
  pValue = test$p.value
  outTab = rbind(outTab, cbind(rt1, gene=i))
}

write.table(outTab,file="all_data.txt",sep="\t",row.names=F,quote=F)


#plot all
data=read.table("all_data.txt",sep="\t",header=T,check.names=F)
data<-na.omit(data)
colnames(data)[2]<-'Group'
data$Group=factor(data$Group, levels=c("control", "UC"))
data$expression<-log2(data$expression)
#data<-data[data$gene != 'T cells CD4 memory activated',]
#data<-data[data$gene != 'B cells memory',]
#data<-data[data$gene != 'T cells CD8',]
#data<-data[data$gene != 'T cells follicular helper',]
#data<-data[data$gene != 'T cells gamma delta',]
#data<-data[data$gene != 'NK cells activated',]
#data<-data[data$gene != 'Dendritic cells resting',]
#data<-data[data$gene != 'Dendritic cells activated',]
#data<-data[data$gene != 'Mast cells resting',]
#data<-data[data$gene != 'Mast cells activated',]
#data<-data[data$gene != 'Neutrophils',]

#第二次删减
#data<-data[data$gene != 'B cells naive',]
#data<-data[data$gene != 'T cells CD4 naive',]
#data<-data[data$gene != 'Macrophages M1',]
#data<-data[data$gene != 'Mast cells resting',]
#data<-data[data$gene != 'Mast cells activated',]

p<-ggboxplot(data, x="gene", y="expression", fill = "Group",
            ylab="Fraction",
            xlab="",
            palette = c("#4DBBD5FF","#DC0000FF"))
p=p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Group),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif", size = 5)+
  theme(text = element_text(size = 12))
p
ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 8, height = 6)


###############################################################################
## 目的：
##   1. 用 rmcorr 重新计算 CCL11/CLDN1 与 22 种免疫细胞的相关性（患者内关联）
##   2. 患者水平聚合 Spearman 作为敏感性分析
##   3. 输出主表、补充表、Figure 6B 热图、Figure 6C/D 散点图
## 输入：
##   uniq.symbol.txt           行=基因，列=样本
##   CIBERSORT-Results.txt     行=样本，列=22 细胞 + P-value/Correlation/RMSE
##   clinical_metadata.txt     含 id、Patient_ID、Disease
##   venn.txt                  两行：CCL11、CLDN1
## 输出：
##   rmcorr_gene_immune.csv                完整 rmcorr 结果（44 行）
##   patient_level_spearman.csv            患者水平聚合 Spearman（44 行）
##   Supplementary_Table_X_correlation.csv 合并后的补充表
##   Figure6B_rmcorr_heatmap.pdf
##   Figure6C_D_rmcorr_scatter.pdf
###############################################################################

rm(list = ls())

## ---------------------------------------------------------------------------
## 0. 路径与依赖
## ---------------------------------------------------------------------------
work_dir <- "D:/LZ-25234/UC/11Immue/05cor"   # 改成你的目录
setwd(work_dir)

if (!requireNamespace("rmcorr", quietly = TRUE)) install.packages("rmcorr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggcorrplot", quietly = TRUE)) install.packages("ggcorrplot")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(rmcorr)
library(ggplot2)
library(ggcorrplot)
library(reshape2)

## ---------------------------------------------------------------------------
## 1. 读表达矩阵，取 CCL11 / CLDN1
## ---------------------------------------------------------------------------
exp <- read.table("uniq.symbol.txt", header = TRUE, row.names = 1,
                  sep = "\t", check.names = FALSE)
exp <- as.matrix(exp)
exp <- t(exp)   # 样本 × 基因

gene_df <- data.frame(
  id    = rownames(exp),
  CCL11 = as.numeric(exp[, "CCL11"]),
  CLDN1 = as.numeric(exp[, "CLDN1"]),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------------
## 2. 读 CIBERSORT，删除统计列
## ---------------------------------------------------------------------------
cibe <- read.table("CIBERSORT-Results.txt", header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE)

drop_cols <- intersect(c("P-value", "Correlation", "RMSE"), colnames(cibe))
if (length(drop_cols) > 0) {
  cibe <- cibe[, setdiff(colnames(cibe), drop_cols), drop = FALSE]
}
cibe$id <- rownames(cibe)
cell_types <- setdiff(colnames(cibe), "id")

cat("CIBERSORT 细胞类型数：", length(cell_types), "\n")
#CIBERSORT 细胞类型数： 22 
cat("CIBERSORT id 示例：", paste(head(cibe$id, 3), collapse = ", "), "\n")

## ---------------------------------------------------------------------------
## 3. 读临床表，拿 Patient_ID，只保留 UC
## ---------------------------------------------------------------------------
clin <- read.table("clinical_metadata.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
clin$id <- as.character(clin$id)

## 如果 clinical_metadata.txt 不在当前目录，改成绝对路径，例如：
## clin <- read.table("D:/LZ-25234/UC/12Subtype/03immue/04cibersort/clinical_metadata.txt",
##                    header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 4. 合并，只保留 UC
## ---------------------------------------------------------------------------
dat <- merge(gene_df, cibe, by = "id")
dat <- merge(dat, clin[, c("id", "Patient_ID", "Disease")], by = "id")
dat <- dat[dat$Disease == "UC", ]
dat$Patient_ID <- factor(dat$Patient_ID)

cat("UC 活检数：", nrow(dat),
    "，患者数：", length(unique(dat$Patient_ID)), "\n")
#UC 活检数： 865 ，患者数： 419

## 检查每个基因/细胞是否有缺失
if (any(is.na(dat$CCL11)) || any(is.na(dat$CLDN1))) {
  warning("CCL11/CLDN1 存在 NA，将在逐对分析中按 complete.cases 剔除")
}

## ---------------------------------------------------------------------------
## 5. rmcorr：2 基因 × 22 细胞
## ---------------------------------------------------------------------------
genes <- c("CCL11", "CLDN1")

rm_list <- list()
for (g in genes) {
  for (ct in cell_types) {
    d_sub <- dat[, c("Patient_ID", g, ct)]
    colnames(d_sub) <- c("Participant", "Measure1", "Measure2")
    d_sub <- d_sub[complete.cases(d_sub), ]
    
    ## rmcorr 要求每个 participant 至少 2 个观测
    keep_pat <- names(which(table(d_sub$Participant) >= 2))
    d_sub <- d_sub[d_sub$Participant %in% keep_pat, ]
    d_sub$Participant <- droplevels(factor(d_sub$Participant))
    
    n_pat <- length(unique(d_sub$Participant))
    n_bio <- nrow(d_sub)
    
    if (n_pat < 3) {
      rm_list[[length(rm_list) + 1]] <- data.frame(
        Gene = g, Cell = ct,
        rmcorr_r = NA, rmcorr_df = NA, rmcorr_p = NA,
        n_patients = n_pat, n_biopsies = n_bio,
        Note = "患者数不足"
      )
      next
    }
    
    fit <- tryCatch(
      rmcorr(participant = Participant, measure1 = Measure1, measure2 = Measure2,
             data = d_sub),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      rm_list[[length(rm_list) + 1]] <- data.frame(
        Gene = g, Cell = ct,
        rmcorr_r = NA, rmcorr_df = NA, rmcorr_p = NA,
        n_patients = n_pat, n_biopsies = n_bio,
        Note = "rmcorr 拟合失败"
      )
    } else {
      rm_list[[length(rm_list) + 1]] <- data.frame(
        Gene = g, Cell = ct,
        rmcorr_r  = as.numeric(fit$r),
        rmcorr_df = as.numeric(fit$df),
        rmcorr_p  = as.numeric(fit$p),
        n_patients = n_pat, n_biopsies = n_bio,
        Note = ""
      )
    }
  }
}

rm_res <- do.call(rbind, rm_list)
rm_res$rmcorr_FDR <- p.adjust(rm_res$rmcorr_p, method = "BH")

## 按基因、细胞排序
rm_res <- rm_res[order(rm_res$Gene, rm_res$Cell), ]
write.csv(rm_res, "rmcorr_gene_immune.csv", row.names = FALSE)
cat("\n===== rmcorr 结果 =====\n")
print(rm_res)

## ---------------------------------------------------------------------------
## 6. 患者水平聚合 Spearman（敏感性分析）
## ---------------------------------------------------------------------------
pat_dat <- aggregate(
  dat[, c("CCL11", "CLDN1", cell_types)],
  by = list(Patient_ID = dat$Patient_ID),
  FUN = mean, na.rm = TRUE
)

cat("\n患者水平聚合纳入患者数：", nrow(pat_dat), "\n")
#患者水平聚合纳入患者数： 419

agg_list <- list()
for (g in genes) {
  for (ct in cell_types) {
    ct_test <- tryCatch(
      cor.test(pat_dat[[g]], pat_dat[[ct]], method = "spearman"),
      error = function(e) NULL
    )
    if (is.null(ct_test)) {
      agg_list[[length(agg_list) + 1]] <- data.frame(
        Gene = g, Cell = ct,
        patient_rho = NA, patient_p = NA,
        n_patients_agg = nrow(pat_dat)
      )
    } else {
      agg_list[[length(agg_list) + 1]] <- data.frame(
        Gene = g, Cell = ct,
        patient_rho = as.numeric(ct_test$estimate),
        patient_p   = as.numeric(ct_test$p.value),
        n_patients_agg = nrow(pat_dat)
      )
    }
  }
}

agg_res <- do.call(rbind, agg_list)
agg_res$patient_FDR <- p.adjust(agg_res$patient_p, method = "BH")
agg_res <- agg_res[order(agg_res$Gene, agg_res$Cell), ]
write.csv(agg_res, "patient_level_spearman.csv", row.names = FALSE)
cat("\n===== 患者水平聚合 Spearman 结果 =====\n")
print(agg_res)

## ---------------------------------------------------------------------------
## 7. 合并成 Supplementary Table X
## ---------------------------------------------------------------------------
supp_tab <- merge(
  rm_res[, c("Gene", "Cell", "rmcorr_r", "rmcorr_p", "rmcorr_FDR",
             "n_patients", "n_biopsies", "Note")],
  agg_res[, c("Gene", "Cell", "patient_rho", "patient_p",
              "patient_FDR", "n_patients_agg")],
  by = c("Gene", "Cell"),
  all = TRUE
)

supp_tab <- supp_tab[order(supp_tab$Gene, supp_tab$Cell), ]
write.csv(supp_tab, "Supplementary_Table_X_correlation.csv", row.names = FALSE)
cat("\n===== Supplementary Table X 已输出 =====\n")
print(supp_tab)

## ---------------------------------------------------------------------------
## 8. Figure 6B：rmcorr 热图
## ---------------------------------------------------------------------------
## 组装 2 × 22 矩阵（即使某些细胞失败也保留 NA）
rm_mat <- dcast(rm_res, Gene ~ Cell, value.var = "rmcorr_r")
rownames(rm_mat) <- rm_mat$Gene
rm_mat <- as.matrix(rm_mat[, -1, drop = FALSE])

p_mat <- dcast(rm_res, Gene ~ Cell, value.var = "rmcorr_FDR")
rownames(p_mat) <- p_mat$Gene
p_mat <- as.matrix(p_mat[, -1, drop = FALSE])

## 对齐
cells_use <- colnames(rm_mat)
rm_mat <- rm_mat[c("CCL11", "CLDN1"), cells_use, drop = FALSE]
p_mat  <- p_mat[c("CCL11", "CLDN1"), cells_use, drop = FALSE]

pdf("Figure6B_rmcorr_heatmap.pdf", width = 6, height = 7)
ggcorrplot(
  t(rm_mat),
  show.legend = TRUE,
  sig.level = 0.05,
  lab = FALSE,
  p.mat = t(p_mat),
  hc.order = FALSE
) + coord_flip() +
  ggtitle("Repeated-measures correlation (rmcorr)\nCCL11/CLDN1 vs immune cells")
dev.off()

## ---------------------------------------------------------------------------
## 9. Figure 6C/D：巨噬细胞 M0/M1/M2 散点图（标注 rmcorr 结果）
## ---------------------------------------------------------------------------
plot_gene_cell <- function(gene, cell, dat, rm_res) {
  d <- dat[, c("Patient_ID", gene, cell)]
  colnames(d) <- c("Patient_ID", "Gene", "Cell")
  d <- d[complete.cases(d), ]
  
  stat <- rm_res[rm_res$Gene == gene & rm_res$Cell == cell, ]
  
  ggplot(d, aes(x = Cell, y = Gene)) +
    geom_point(alpha = 0.4, color = "#4C72B0", size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    labs(
      x = paste0(cell, " infiltration"),
      y = paste0(gene, " expression"),
      title = paste0(gene, " vs ", cell),
      subtitle = sprintf(
        "rmcorr r = %.2f, P = %.2g (n = %d patients, %d biopsies)",
        stat$rmcorr_r, stat$rmcorr_p, stat$n_patients, stat$n_biopsies
      )
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

pdf("Figure6C_D_rmcorr_scatter.pdf", width = 10, height = 3.5)
print(plot_gene_cell("CCL11", "Macrophages M0", dat, rm_res))
print(plot_gene_cell("CCL11", "Macrophages M1", dat, rm_res))
print(plot_gene_cell("CCL11", "Macrophages M2", dat, rm_res))
print(plot_gene_cell("CLDN1", "Macrophages M0", dat, rm_res))
print(plot_gene_cell("CLDN1", "Macrophages M1", dat, rm_res))
print(plot_gene_cell("CLDN1", "Macrophages M2", dat, rm_res))
dev.off()

cat("\n===== 全部完成 =====\n")
cat("输出文件：\n")
cat("  rmcorr_gene_immune.csv\n")
cat("  patient_level_spearman.csv\n")
cat("  Supplementary_Table_X_correlation.csv\n")
cat("  Figure6B_rmcorr_heatmap.pdf\n")
cat("  Figure6C_D_rmcorr_scatter.pdf\n")




rm(list=ls())

path <- 'D:\\LZ-25234\\UC\\11Immue\\05cor'
setwd(path)

library(ggplot2)

## ---------------------------------------------------------------------------
## 1. 读表达矩阵与 CIBERSORT
## ---------------------------------------------------------------------------
exp <- read.table('uniq.symbol.txt', header = TRUE, row.names = 1, sep = '\t')
cibe <- read.table('CIBERSORT-Results.txt', header = TRUE, row.names = 1, sep = '\t')

exp <- t(exp)
exp <- data.frame(exp)

## 如果 CIBERSORT 列名是点号（Macrophages.M0），保持原样；
## 如果列名是空格，请把下面 gsub 改成 gsub("\\.", " ", ...)
# colnames(cibe) <- gsub("\\.", " ", colnames(cibe))

## ---------------------------------------------------------------------------
## 2. 读 rmcorr 结果
## ---------------------------------------------------------------------------
rm_res <- read.csv("rmcorr_gene_immune.csv",
                   check.names = FALSE,
                   stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 3. 定义画图函数（视觉风格与原 scatte.R 一致）
## ---------------------------------------------------------------------------
plot_scatter_rmcorr <- function(gene, cell_rm, cell_cibe, rm_res) {
  df <- data.frame(
    Immune = cibe[[cell_cibe]],
    Gene   = exp[[gene]]
  )
  df <- df[complete.cases(df), ]
  
  stat <- rm_res[rm_res$Gene == gene & rm_res$Cell == cell_rm, ]
  if (nrow(stat) == 0) stop("rmcorr 结果中找不到：", gene, " - ", cell_rm)
  
  label_text <- sprintf("rmcorr r = %.2f, P = %.2e", stat$rmcorr_r, stat$rmcorr_p)
  
  ## 动态放置统计标注，避免遮挡数据
  x_min <- min(df$Immune, na.rm = TRUE)
  x_max <- max(df$Immune, na.rm = TRUE)
  y_min <- min(df$Gene, na.rm = TRUE)
  y_max <- max(df$Gene, na.rm = TRUE)
  label_x <- x_min + 0.05 * (x_max - x_min)
  label_y <- y_max - 0.15 * (y_max - y_min)
  
  p <- ggplot(df, aes(x = Immune, y = Gene)) +
    geom_point(color = "#2C7BB6", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "#D7191C", se = TRUE, linetype = "dashed") +
    annotate("text", x = label_x, y = label_y, label = label_text,
             size = 4.5, hjust = 0, color = "black") +
    theme_minimal() +
    labs(
      title = cell_rm,
      x = "Immune cell infiltration level",
      y = paste0(gene, " expression level")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  return(p)
}

## ---------------------------------------------------------------------------
## 4. 批量生成 6 张图（CCL11/CLDN1 × M0/M1/M2）
## ---------------------------------------------------------------------------
genes <- c("CCL11", "CLDN1")
cells <- c("Macrophages M0", "Macrophages M1", "Macrophages M2")

for (g in genes) {
  for (cl in cells) {
    cl_cibe <- gsub(" ", ".", cl)  # rmcorr 里是空格，CIBERSORT 列名是点号
    p <- plot_scatter_rmcorr(g, cl, cl_cibe, rm_res)
    
    file_name <- paste0(g, "_", gsub(" ", "_", cl), "_rmcorr.pdf")
    pdf(file_name, height = 3.5, width = 3.5)
    print(p)
    dev.off()
    cat("Saved:", file_name, "\n")
  }
}
