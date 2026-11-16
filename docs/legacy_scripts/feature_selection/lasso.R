

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
