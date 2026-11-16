

rm(list=ls())

#path
path<-'E:\\LZ-25007\\colitis\\06verify\\GSE87466\\01model'
setwd(path)

#data
rt<-read.table('normalize.txt',header = T,row.names=1,sep='\t')
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')

rt<-rt[hub$V1,]
rt<-t(rt)
rt<-data.frame(rt)

rt$Sample<-rownames(rt)

rt<-merge(group,rt)
rownames(rt)<-rt[,1]
rt<-rt[,-1]
colnames(rt)[1]<-'Group'
table(rt$Group)
#control      UC 
#     21      87 


#训练最终模型
rt1<-rt
dim(rt1)
#[1] 108   3

head(rt1)

#处理分组
rt1[which(rt1$Group == 'UC'),]$Group<-1
rt1[which(rt1$Group == 'control'),]$Group<-0
rt1$Group<-as.numeric(rt1$Group)

#ROC
f <- glm(Group~CCL11+MMP1, data = rt1, family = binomial())
summary(f)
glm<-data.frame(f$coefficients)
#pred <- predict(f, newdata = rt1, type = "response")


library(ROCR)
prob <- predict(f, newdata = rt1, type = "response")
pred <- prediction(prob, rt1$Group)
performance(pred, "auc")@y.values[[1]]  #[1] 0.9846743
#plot(performance(pred, "tpr", "fpr"), colorize = F, lwd = 3, main = "ROC Curves")
rt1<-cbind(rt1,prob)

#plot
library(timeROC)
library(pROC)
roc1<- roc(rt1$Group, rt1$prob)
roc2<- roc(rt1$Group, rt1$CCL11)
roc3<- roc(rt1$Group, rt1$MMP1)

#dev.new()
plot(roc1,col="red")
plot(roc2,add=TRUE,  col="green")
plot(roc3, add=TRUE, col="blue")
legend('bottomright',
       c(paste0('GSE87466: ', round(roc1$auc, 3)),
         paste0('CCL11: ', round(roc2$auc, 3)),
         paste0('MMP1: ', round(roc3$auc, 3))
       ),
       col = c("red", "green", "blue"),
       lwd = 2, bty = 'n')

#保存模型
tiff(file="model.tiff", width =4, height =4, units ="in", bg="white",
     res= 350)
#pdf(file="model.pdf", width =4, height =4)
plot(roc1,col="red")
plot(roc2,add=TRUE,  col="green")
plot(roc3, add=TRUE, col="blue")
legend('bottomright',
       c(paste0('GSE87466: ', round(roc1$auc, 3)),
         paste0('CCL11: ', round(roc2$auc, 3)),
         paste0('MMP1: ', round(roc3$auc, 3))
       ),
       col = c("red", "green", "blue"),
       lwd = 2, bty = 'n')
dev.off()

