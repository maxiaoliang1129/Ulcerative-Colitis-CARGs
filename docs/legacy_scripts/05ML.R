

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\05ML\\00lasso'
setwd(path)

#data
diff_expr<-read.table('diffGeneExp_p0.05.txt',header=T,row.names=1,sep='\t')
diff_expr<-t(diff_expr)
dim(diff_expr)
#[1] 1320  779


#分组文件
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Group')

diff_expr<-diff_expr[group$Sample,]
dim(diff_expr);dim(group)
#[1] 1320  779
#[1] 1320    2

#这里保留交集
gene<-read.table('diff_intergene.txt',header=F,sep='\t')
gene<-gene$V1

#这里把表达矩阵的基因全部替换成大写
diff_expr<-diff_expr[,gene]
dim(diff_expr)
#[1] 1320   35

table(group$Group)
#control  UC
#455     865

#glm
x<-diff_expr
group[which(group$Group == 'UC'),]$Group<-1
group[which(group$Group == 'control'),]$Group<-0
y<-factor(group$Group)

df2<-cbind(group$Group,diff_expr)
colnames(df2)[1]<-'Group'
df2<-data.frame(df2)
df2[,2:ncol(df2)]<-as.numeric(unlist(df2[,2:ncol(df2)]))
df2$Group<-factor(df2$Group)
dim(df2)
#[1] 1320   36

#共线性诊断
XX<-cor(df2[-1])
kappa(XX,exact = TRUE)
#[1] 260.5738
#存在共线性 lasso回归解决共线性问题

#-----------------------------

library(glmnet)

#使用整套做训练集
x <- as.matrix(df2[,-1])
y <- df2$Group
y<-as.numeric(y)

#Lasso回归
f1<-glmnet(x,y,family="binomial",nlambda=100,alpha=1)
print(f1)

n<-min(f1$lambda)
n
#[1] 0.0001843275

#输出图形
plot(f1, xvar="lambda", label=TRUE)

#save
pdf("lambda.pdf",width = 4, height =4)
plot(f1, xvar="lambda", label=TRUE)
dev.off()


#交叉验证
cvfit<-cv.glmnet(x,y)
plot(cvfit)

#save
pdf("cvfit.pdf",width = 4, height =4)
plot(cvfit)
dev.off()


#选择最佳lambada值构建最佳模型
#cvfit$lambda.min#求出最小值一个标准误的λ值
cvfit$lambda.1se#求出最小值
#[1] 0.001427181
#[1] 0.006939814

#带入模型
l.coef1<-coef(cvfit$glmnet.fit,s=cvfit$lambda.1se,exact = F)
#l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact = F)
l.coef1
#l.coef2


#save 最后结果
coef <- coef(f1,s = cvfit$lambda.1se)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene = row.names(coef)[index]
geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
dim(geneCoef) #[1] 25  2

write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names = F)




rm(list=ls())

#library(xlsx)
library(randomForest)
library(varSelRF)
library(pROC)

library(rfPermute)
#library(rfUtilities)

#设置种子数
set.seed(123) #固定结果

#path
path<-'D:\\LZ-25234\\UC\\05ML\\01RF'
setwd(path)

diff_expr<-read.table('diffGeneExp_p0.05.txt',header=T,row.names=1,sep='\t')
diff_expr<-t(diff_expr)
dim(diff_expr)
#[1] 1320  779


#分组文件
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Group')

diff_expr<-diff_expr[group$Sample,]
dim(diff_expr);dim(group)
#[1] 1320  779
#[1] 1320    2

#这里保留交集
gene<-read.table('diff_intergene.txt',header=F,sep='\t')
gene<-gene$V1

#这里把表达矩阵的基因全部替换成大写
rowsid<-colnames(diff_expr)
rowsid<-toupper(rowsid)
colnames(diff_expr)<-rowsid

diff_expr<-diff_expr[,gene]
dim(diff_expr)
#[1] 1320   35


#rf
x<-diff_expr
y<-factor(group$Group)

df2<-cbind(group$Group,diff_expr)
colnames(df2)[1]<-'Group'

#设置种子
set.seed(123)
#随机森林树
rf<-randomForest(x,y,data=df2,importance=TRUE,proximity=TRUE,ntree=800)
rf
plot(rf)

plot(margin(rf, df2$Group), main = '观测值被判断正确的概率图')

which.min(rf$err.rate[,1])
#[1] 118

#save 随机森林树
pdf('forest.pdf',width=4,height=4)
plot(rf,main="Random forest",lwd=2)
dev.off()

#获取最小决策树后 再次构建模型
set.seed(123)
rf2<-randomForest(x,y,data=df2,importance=TRUE,proximity=TRUE,ntree=118)
rf2

getTree(rf2,1)


##关键特征基因识别
#查看表示每个变量重要性的得分
#summary(rf2)
importance_gene<-rf2$importance
head(importance_gene)

'''
#或者使用函数importance()
importance_gene<-data.frame(importance(rf2))
head(importance_gene)
'''

#可以根据某种重要性的高低排个序,例如根据“MeanDecreaseGini”指标
importance_gene<-importance_gene[order(importance_gene[,4], decreasing = TRUE), ]
head(importance_gene)
 
#输出表格
write.table(importance_gene, 'importance_gene.txt', sep = '\t', col.names = NA, quote = FALSE)

#查看“MeanDecreaseGini”指标>15的基因
nums<-data.frame(importance_gene)
nums<-nums[which(nums$MeanDecreaseGini > 15),]
dim(nums)
#[1] 19  4

#作图展示 top14 重要的 OTUs
varImpPlot(rf2, 
     n.var = min(35, nrow(rf2$importance)), 
     main = 'Total 35 - variable importance')
#保存图片
pdf(file="VIP.pdf", width =6, height =8)
varImpPlot(rf2, 
           n.var = min(35, nrow(rf2$importance)), 
           main = 'Total 35 - variable importance')
dev.off()


#交叉验证寻找合适数量的变量
#5次重复十折交叉验证
set.seed(123)
df2<-as.matrix(df2)
otu_train.cv<-replicate(5,rfcv(x, y, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv
 
#提取验证结果绘图
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$gene <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'gene')
otu_train.cv$gene <- as.numeric(as.character(otu_train.cv$gene))
 
#拟合线图
library(ggplot2)
library(splines)
#用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线
 
p <- ggplot(otu_train.cv, aes(gene, value)) +
    geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of genes', y = 'Cross-validation error')
 
p

pdf("cross.pdf", onefile=F,width =4, height =4)
p
dev.off()


#挑选top3个打分的gene
#选择top3重要的genes MeanDecreaseGini>15
otu_select<-rownames(importance_gene)[1:19]
otu_select

#导出基因
write.table(data.frame(otu_select),'rf.genes.txt',sep='\t',
     col.names=F,row.names=F,quote=F)


##===================================
'''
##寻找重要变量
#重要性评分
importance_rf<-importance(rf2,scale = TRUE)
importance_otu.scale<-data.frame(importance_rf,check.names = FALSE)
importance_otu.scale

dim(importance_rf)
#[1] 14  4

#根据评分排序
#对预测变量按重要性得分排个序，例如根据“MeanDecreaseGini”
importance_otu.scale<-importance_otu.scale[order(importance_otu.scale$MeanDecreaseGini, decreasing = TRUE), ]

#save 保存评分结果
write.table(importance_rf,'importance_rf.txt',row.names = T,col.names = T,
            sep='\t',quote = F)

#绘图
#简单地作图展示预测变量的MeanDecreaseGini值
library(ggplot2)
 
importance_otu.scale$gene<-rownames(importance_otu.scale)
importance_otu.scale$gene<-factor(importance_otu.scale$gene, levels = importance_otu.scale$gene)
 
p <- ggplot(importance_otu.scale, aes(gene, MeanDecreaseGini)) +
geom_col(width = 0.5, fill = '#FFC068', color = NA) +
labs(title = NULL, x = NULL, y = 'MeanDecreaseGini', fill = NULL) +
theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_y_continuous(expand = c(0, 0), limit = c(0, 16))
 
p
 
#右上角备注模型的已知解释率
p <- p +
annotate('text', label = 'Plant Age', x = 9, y = 15, size = 4) +
annotate('text', label = sprintf('italic(R^2) == %.2f', 96.14), x = 9, y = 13, size = 3, parse = TRUE)
 
p

'''

##===================================

'''
#VIP得分展示
#varImpPlot(rf,main = "variable importance")
varImpPlot(rf3,sort=TRUE,n.var=nrow(head(rf3$importance,10)),
           main = "variable importance")

#save plot
pdf("geneImporttance.pdf",width=6.2,height=6)
varImpPlot(rf3,sort=TRUE,n.var=nrow(head(rf3$importance,10)),
           main = "variable importance")
dev.off()

tiff("geneImporttance.tiff",width=6.2,height=6,units = "in", bg = "white",
     res = 500)
varImpPlot(rf3,sort=TRUE,n.var=nrow(head(rf3$importance,10)),
           main = "variable importance")
dev.off()


#筛选特征基因
rfGenes<-importance_rf[order(importance_rf[,"MeanDecreaseGini"],decreasing=TRUE),]

'''
#VIP>2
rfGenes<-as.data.frame(rfGenes)
rfGenes<-rownames(rfGenes[rfGenes$MeanDecreaseGini>2,])
'''
#VIP评分最高的10个基因
rfGenes<-rownames(rfGenes[1:10,])

#show
rfGenes
#[1] "PNRC1"  "HS3ST1" "ZFP36"  "TES"    "BCAN"   "IER3"   "CDKN1A" "IGFBP1" "CHST2"  "CP"

#save
write.table(rfGenes,'rfGenes.txt',sep='\t',quote=F,col.names=F,row.names=F)

'''



#输出表达量  tpm矩阵
sigExp<-expr[rfGenes,]

#save
write.table(sigExp,'rfGenes_expr.txt',row.names=T,col.names=T,sep='\t',quote=F)
'''

-----------------------------------


#预测：用测试集
test_data<-read.table('testdata1.txt',header = T,row.names = 1,sep='\t')
testdata<-test_data
colnames(testdata)[14]<-'NME1-NME2'
#testdata<-testdata[,-1]
#testdata<-df2

iris_pred<-predict(rf3, newdata=testdata)

testdata<-data.frame(testdata)
table(iris_pred);table(testdata$Group)
#iris_pred
#Normal  Tumor 
#     5    580 

#Normal  Tumor 
#    59    526 

'''
#模型评价AUC
iris_pred<-data.frame(iris_pred)
roc<-multiclass.roc(as.ordered(test_data$Group) ,as.ordered(iris_pred$iris_pred))
roc

roc<-roc(as.ordered(testdata$Group) ,as.ordered(iris_pred$iris_pred))
plot(roc,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2),
     grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
'''




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\05ML\\03XGBoost'
setwd(path)

#BiocManager::install('xgboost')
library(xgboost)

#data
diff_expr<-read.table('diffGeneExp_p0.05.txt',header=T,row.names=1,sep='\t')
diff_expr<-t(diff_expr)
dim(diff_expr)
#[1]  124 1024

#分组文件
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Group')

diff_expr<-diff_expr[group$Sample,]
dim(diff_expr);dim(group)

#这里保留交集
gene<-read.table('diff_intergene.txt',header=F,sep='\t')
gene<-gene$V1

diff_expr<-diff_expr[,gene]
dim(diff_expr)

table(group$Group)

x<-diff_expr
group[which(group$Group == 'UC'),]$Group<-1
group[which(group$Group == 'control'),]$Group<-0
y<-as.numeric(group$Group)

dim(x);length(y) #y的数量必须和x的行数一致
#[1] 1320   35
#[1] 1320

#XGBoost
set.seed(1234)
dtrain <- xgb.DMatrix(data = as.matrix(x), label = y)

##2.设置参数，二分类任务为例
params <- list(
  objective = "binary:logistic", 
  eval_metric = "auc", 
  max_depth = 6, 
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

##3.训练模型


#自动选择最佳树
cv <- xgb.cv(params = params, data = dtrain, nrounds = 500, nfold = 10,
             early_stopping_rounds = 10, verbose = 0, metrics = "auc")
best_nrounds <- cv$best_iteration
best_nrounds
#[1] 94
nrounds<-best_nrounds

#nrounds <- 100
model <- xgb.train(params = params, data = dtrain, nrounds = best_nrounds, verbose = 1)

#4.提取特征重要性
importance_matrix <- xgb.importance(feature_names = colnames(x), model = model)

#5.查看重要性表的前几行
head(importance_matrix)

write.table(importance_matrix,'importance_matrix.txt',row.names = F,col.names = T,sep = '\t',quote = F)

#重要性指标含义：Gain（提升模型的贡献度）、Cover（平均覆盖率）、Frequency（使用频率）
#选择最高Gain值的基因通常代表更关键的生物学标记。

#6.选择关键特征基因（例如选Top 10个）
top_genes <- importance_matrix$Feature[1:10]
print(top_genes)
# [1] "S100A9"  "S100A12" "IL1R2"   "GSTO1"   "UQCRH"   "PRF1"    "TBX21"   "RAB32"   "LBH"     "S1PR5"  

#保存结果
write.table(data.frame(top_genes),'xgboost.gene.txt',row.names = F,col.names = F,sep = '\t',quote = F)

pdf(file = "xgboost.pdf",height = 4,width = 4.5)
xgb.plot.importance(importance_matrix[1:10,])
dev.off()
#或者
#barplot(importance_matrix[,1])




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\05ML\\04venn'
setwd(path)

#data
rf<-read.table('rf.genes.txt',sep='\t')
lasso<-read.table('lasso.genes.txt',sep='\t')
svm<-read.table('xgboost.gene.txt',sep='\t')

a<-rf$V1
b<-lasso$V1
c<-svm$V1

length(a);length(b);length(c)
#[1] 19
#[1] 25
#[1] 10

#绘图
library(VennDiagram)
gene_list<-list('Random Forest'=a,'LASSO'=b,'XGBoost'=c)
#save upset
v1 <- venn.diagram(x = gene_list, 
                   filename = NULL, # 直接给一个名称会自动保存文件到本地
                   scaled = F, # 根据比例显示大小
                   alpha=c(0.3, 0.3,0.3), # 设置每个区块的透明度
                   lwd = 1, # 圆圈线条的粗细：1 2 3 4 5 6
                   lty = 'blank', # 圆圈线条的类型：1为实线，2为虚线，blank为无线条
                   col = c("black","red",'blue'), # 圆圈线条颜色
                   fill = c("skyblue", "springgreen",'lightpink'),  # 圆圈颜色
                   cex = 1, # 数字大小
                   fontface = "bold", # 加粗
                   fonrfamily = "sans", # 数字字体
                   cat.cex = 1,  # 标签字体大小
                   cat.col = "black",  # 标签字体色
                   cat.fontface = "bold",  # 加粗
                   cat.default.pos = "outer",  # 标签内外位置, 在圆圈内还是圆圈外，outer 内 text 外
                   cat.pos = c(0, 0, 0),  # 标签旋转位置，用圆的度数
                   cat.dist = c(0.05,0.05,0.05),  # 标签离圆圈位置,离圆的距离，如果标签与圆圈重叠，可以调整这个参数
                   cat.fontfamily = "sans",  # 标签字体
)
cowplot::plot_grid(v1)

pdf("Venn.pdf", onefile=F, width =4, height =4)
cowplot::plot_grid(v1)
dev.off()

#接下来把至少出现在两种过滤方法中的基因保留下来
d1<-intersect(a,b)
d1
#[1] "S100A9"  "S100A12" "IL1R2"   "GSTO1"   "UQCRH"  


d2<-intersect(a,c)
d2
#[1] "UQCRH" "LBH"  

d3<-intersect(b,c)
d3
#[1] "UQCRH" "SYNE2"

d<-c(d1,d2,d3)
d<-unique(d)
d;length(d)
#[1] "S100A9"  "S100A12" "IL1R2"   "GSTO1"   "UQCRH"   "LBH"     "SYNE2"  
#[1] 7

#保存结果
write.table(data.frame(d),'venn.txt',row.names = F,col.names = F,sep = '\t',quote = F)


##plot2
library(UpSetR)

p<-list('Random Forest'=a,'LASSO'=b,'XGBoost'=c)

#save plot
pdf("upset.pdf", onefile=F, width =5, height =4)
upset(fromList(p), sets.bar.color=c("red", "blue","yellow"),
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2))
dev.off()

#交集分析
x<-intersect(a,b)
d<-intersect(x,c)
length(d) #[1] 9
d
write.table(data.frame(d),'venn.txt',row.names = F,col.names = F,sep = '\t',quote = F)





