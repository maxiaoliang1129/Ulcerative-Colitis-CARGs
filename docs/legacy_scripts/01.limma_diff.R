

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\01diff'
setwd(path)

library(limma)
library(impute)
library(edgeR)
library(DESeq2)

#data
exp<-read.table("normalize.txt", sep = "\t", header = T)
exp<-as.matrix(exp)
dim(exp)
#[1] 39374  1320

dimnames<-list(rownames(exp),colnames(exp))
exp<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt<-exp
dim(rt)
#[1] 39374  1320

####差异分析
#method1
#------------------------------------
#分组信息
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Type')
dim(group) #[1] 1320    2

rt<-grouprt<-rt[,group$Sample]

rownames(group)<-group$Sample
sample<-group$Sample
group<-data.frame(group$Type)
rownames(group)<-sample
colnames(group)[1]<-'Group'
group$Group<-as.factor(group$Group)
group<-group$Group
design<-model.matrix(~0+group)

rownames(design)<-colnames(rt)
head(design)

#无需进行voom转换
#线性模型拟合
exprSet<-rt
fit<-lmFit(exprSet, design)

contrast.matrix <- makeContrasts(groupUC - groupcontrol, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, adjust.method="BH", number=Inf)
dim(deg)  #[1] 39374     6

#Filtering thresholds
adjP=0.05
fc=0.585

output<-deg
#result1-fc=0.584
outdiff<-output[(output$adj.P.Val < adjP & (output$logFC>fc | output$logFC<(-fc))),]
#outdiff<-output[(output$P.Value < adjP),]
dim(outdiff)
#[1] 779    6


#基因上下调
a<-outdiff
a$Group<-'stable'
a[which(a$logFC > 0),'stable']<-'up'
a[which(a$logFC < 0),'stable']<-'down'
table(a$stable)
#down   up 
# 345  434

head(outdiff,10)

write.table(outdiff,'outdiff_P0.05.txt',row.names = T,sep = '\t',quote = F)

#差异基因表达谱
Normalize_counts<-read.table("normalize.txt", sep = "\t", header = T)
outdiff<-read.table('outdiff_P0.05.txt',header = T,sep = '\t')
data<-Normalize_counts
heatmap2=data[rownames(outdiff),]
dim(heatmap2) #[1]  779 1320

#save
write.table(heatmap2,file="diffGeneExp_p0.05.txt",sep="\t",col.names=T,row.names=T,quote=F)


####火山图绘制
outTab<-output
head(outTab)

res<-outTab
logFC<-as.numeric(res$logFC)
FDR<-res$adj.P.Val


#p<0.05 logfc>1
pdf(file="Volcano.pdf", width = 4.5, height = 4.5)
plot(logFC, -1*log10(FDR), 
     col = ifelse(FDR <= 0.05 & logFC > 0.585, "red", 
                  ifelse(FDR <= 0.05 & logFC < -0.585, "blue", "grey")),
     xlab = "logFC", 
     ylab = "-log10(adj.P.Val)", 
     main = "UC'Volcano plot\n(|logFC|>1 & adj.P.value<0.05)", 
     pch = 20,  # 使用实心点，比"."更清晰
     cex = 1.2, # 适当增大点的大小
     xlim = c(-3, 3))

# 添加参考线
abline(v = c(-0.585, 0.585), lty = 2, col = "black")  # logFC阈值线
abline(h = -log10(0.05), lty = 2, col = "black")  # FDR阈值线

# 添加图例
legend("topleft", 
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "grey"), 
       pch = 20,
       bty = "n")
dev.off()




rm(list=ls())

path <- 'D:\\LZ-25234\\UC\\01diff\\敏感性分析'
setwd(path)

library(limma)
library(data.table)


#############################
# 一、主分析：duplicateCorrelation 校正患者内相关性
#    （用 5000 基因子集估计 consensus，再用于全量 lmFit）
#############################

# ---------- 1. 表达矩阵 ----------
exp <- read.table("normalize.txt", sep = "\t", header = T)
exp <- as.matrix(exp)
dimnames <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt  <- exp
dim(rt)   # 39374 1320

# ---------- 2. 分组信息（Sample / Type / Patient） ----------
group_all <- fread('sample_patient_group.txt')
group_all <- as.data.frame(group_all)
stopifnot(all(group_all$Sample %in% colnames(rt)))

# 按 group_all  
rt <- rt[, group_all$Sample]

# ---------- 3. design ----------
group  <- factor(group_all$Type, levels = c("control", "UC"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)          # "control" "UC"
rownames(design) <- colnames(rt)
block <- factor(group_all$Patient)
stopifnot(!anyNA(block), length(block) == ncol(rt))

# ---------- 4. 用 5000 基因子集估计 consensus correlation ----------
set.seed(20240101)
idx <- sample(seq_len(nrow(rt)), 5000)
cat("Estimating consensus correlation on 5000 genes...\n") 
dupcor_sub <- duplicateCorrelation(rt[idx, ], design, block = block)
consensus  <- dupcor_sub$consensus
cat("Consensus correlation =", consensus, "\n")
#Consensus correlation = 0.14215

# 若为负值，保守设为 0
if (consensus < 0) {
  cat("Consensus correlation is negative; setting to 0.\n")
  consensus <- 0
}

write.table(data.frame(consensus_correlation = consensus),
            "dupcor_consensus.txt", row.names = FALSE, quote = FALSE)

# ---------- 5. 用该 consensus 跑全量 lmFit ----------
cat("Fitting limma on full expression matrix...\n")
fit  <- lmFit(rt, design, block = block, correlation = consensus)
contrast.matrix <- makeContrasts(UC - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, number = Inf, adjust.method = "BH")
dim(deg)
#[1] 39374     6

# ---------- 6. 阈值筛选 ----------
adjP <- 0.05
fc   <- 0.585

output  <- deg
outdiff <- output[(output$adj.P.Val < adjP &
                     (output$logFC > fc | output$logFC < -fc)), ]
dim(outdiff)
#[1] 768   6

a <- outdiff
a$Group <- 'stable'
a[which(a$logFC > 0), 'Group'] <- 'up'
a[which(a$logFC < 0), 'Group'] <- 'down'
table(a$Group)
#down   up 
#343  425 

write.table(outdiff, 'outdiff_P0.05_dupcor.txt',
            row.names = TRUE, sep = '\t', quote = FALSE)

# ---------- 7. 差异基因表达谱 ----------
Normalize_counts <- read.table("normalize.txt", sep = "\t", header = T)
data <- Normalize_counts
heatmap2 <- data[rownames(outdiff), ]
dim(heatmap2)
write.table(heatmap2, file = "diffGeneExp_p0.05_dupcor.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

# ---------- 8. 火山图 ----------
res   <- output
logFC <- as.numeric(res$logFC)
FDR   <- res$adj.P.Val

pdf(file = "Volcano_dupcor.pdf", width = 4.5, height = 4.5)
plot(logFC, -1*log10(FDR),
     col = ifelse(FDR <= 0.05 & logFC >  fc, "red",
                  ifelse(FDR <= 0.05 & logFC < -fc, "blue", "grey")),
     xlab = "logFC",
     ylab = "-log10(adj.P.Val)",
     main = paste0("UC vs control (dupcor)\n(|logFC|>", fc, " & adj.P<0.05)"),
     pch = 20, cex = 1.2, xlim = c(-3, 3))
abline(v = c(-fc, fc), lty = 2, col = "black")
abline(h = -log10(0.05), lty = 2, col = "black")
legend("topleft",
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "grey"), pch = 20, bty = "n")
dev.off()


#############################
# 二、敏感性分析：每位患者一个活检（随机版）
#############################
rm(list=ls())
setwd('D:\\LZ-25234\\UC\\01diff\\敏感性分析')

library(limma)
library(data.table)

# ---------- 1. 表达矩阵 ----------
exp <- read.table("normalize.txt", sep = "\t", header = T)
exp <- as.matrix(exp)
dimnames <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt  <- exp
dim(rt)
#[1] 39374  1320

# ---------- 2. 分组信息（每人一个活检，随机版） ----------
group <- fread('group_one_per_patient.txt')
group <- as.data.frame(group)
colnames(group) <- c('Sample', 'Type')
dim(group)                 # 660 2
table(group$Type)          # control 241 / UC 419
#control      UC 
# 241     419 

sum(!group$Sample %in% colnames(rt))   # 必须 0
#[1] 0
stopifnot(all(group$Sample %in% colnames(rt)))

rt <- rt[, group$Sample]
dim(rt)                    # 39374 660
#[1] 39374   660

# ---------- 3. design ----------
rownames(group) <- group$Sample
group$Type <- factor(group$Type, levels = c("control", "UC"))

design <- model.matrix(~ 0 + group$Type)
colnames(design) <- levels(group$Type)   # "control" "UC"
rownames(design) <- colnames(rt)

# ---------- 4. limma 拟合（无重复，不加 block） ----------
fit <- lmFit(rt, design)
contrast.matrix <- makeContrasts(UC - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg  <- topTable(fit2, number = Inf, adjust.method = "BH")
dim(deg)
#[1] 39374     6

# ---------- 5. 阈值筛选 ----------
adjP <- 0.05
fc   <- 0.585

output  <- deg
outdiff <- output[(output$adj.P.Val < adjP &
                     (output$logFC > fc | output$logFC < -fc)), ]
dim(outdiff)

a <- outdiff
a$Group <- 'stable'
a[which(a$logFC > 0), 'Group'] <- 'up'
a[which(a$logFC < 0), 'Group'] <- 'down'
table(a$Group)
#down   up 
#448  557 

write.table(outdiff, 'outdiff_P0.05_onePerPatient.txt',
            row.names = TRUE, sep = '\t', quote = FALSE)

# ---------- 6. 差异基因表达谱 ----------
Normalize_counts <- read.table("normalize.txt", sep = "\t", header = T)
data <- Normalize_counts
heatmap2 <- data[rownames(outdiff), ]
dim(heatmap2) #[1] 1005 1320
write.table(heatmap2, file = "diffGeneExp_p0.05_onePerPatient.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

# ---------- 7. 火山图 ----------
res   <- output
logFC <- as.numeric(res$logFC)
FDR   <- res$adj.P.Val

pdf(file = "Volcano_onePerPatient.pdf", width = 4.5, height = 4.5)
plot(logFC, -1*log10(FDR),
     col = ifelse(FDR <= 0.05 & logFC >  fc, "red",
                  ifelse(FDR <= 0.05 & logFC < -fc, "blue", "grey")),
     xlab = "logFC",
     ylab = "-log10(adj.P.Val)",
     main = paste0("UC vs control (one biopsy/patient)\n(|logFC|>", fc, " & adj.P<0.05)"),
     pch = 20, cex = 1.2, xlim = c(-3, 3))
abline(v = c(-fc, fc), lty = 2, col = "black")
abline(h = -log10(0.05), lty = 2, col = "black")
legend("topleft",
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "grey"), pch = 20, bty = "n")
dev.off()


#############################
# 三、三套结果对比
#############################
deg_orig <- read.table('D:/LZ-25234/UC/01diff/outdiff_P0.05.txt',
                       header = TRUE, sep = '\t', row.names = 1)
deg_dup  <- read.table('outdiff_P0.05_dupcor.txt',
                       header = TRUE, sep = '\t', row.names = 1)
deg_sens <- read.table('outdiff_P0.05_onePerPatient.txt',
                       header = TRUE, sep = '\t', row.names = 1)

c(original      = nrow(deg_orig),
  dupcor        = nrow(deg_dup),
  onePerPatient = nrow(deg_sens))
#     original        dupcor onePerPatient 
#      779           768          1005 

# 原始 vs 每人一样本
common <- intersect(rownames(deg_orig), rownames(deg_sens))
length(common) #[1] 684
mean(sign(deg_orig[common,'logFC']) == sign(deg_sens[common,'logFC']))
#[1] 1
cor(deg_orig[common,'logFC'], deg_sens[common,'logFC'], method = 'spearman')
#[1] 0.9555108

# 校正版 vs 每人一样本
common2 <- intersect(rownames(deg_dup), rownames(deg_sens))
length(common2) #[1] 674
mean(sign(deg_dup[common2,'logFC']) == sign(deg_sens[common2,'logFC'])) #[1] 1
cor(deg_dup[common2,'logFC'], deg_sens[common2,'logFC'], method = 'spearman')
#[1] 0.9561887

# 原始 vs 校正版
common3 <- intersect(rownames(deg_orig), rownames(deg_dup))
length(common3) #[1] 757
mean(sign(deg_orig[common3,'logFC']) == sign(deg_dup[common3,'logFC'])) #[1] 1
cor(deg_orig[common3,'logFC'], deg_dup[common3,'logFC'], method = 'spearman')
#[1] 0.9984491




