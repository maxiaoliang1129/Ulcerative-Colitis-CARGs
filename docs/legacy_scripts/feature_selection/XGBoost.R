

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

