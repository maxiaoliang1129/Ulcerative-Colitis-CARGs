

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\06Model\\04verify\\GSE128682\\02boxplot'
setwd(path)

#data
rt<-read.table('normalize.txt',header = T,sep='\t',row.names=1)
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Group')

rt<-rt[hub$V1,]

rt<-t(rt)
rt<-rt[group$Sample,]
data = data.frame()
rows = na.omit(match(rownames(rt), group$Sample))
table(group$Group)
#control      UC  
#     16      14 
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[,i],gene=i,Group=as.vector(group[rows,2])))
  a[,1] = as.numeric(a[,1])
  test=wilcox.test(expression ~ Group, data = a)
  pValue = test$p.value
  data=rbind(data,a)
}

#save result
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

#plot
library(ggplot2)
library(ggpubr)

data=read.table("data.txt",sep="\t",header=T,check.names=F)

data$Group = factor(data$Group, levels = c("control", "UC"))


##箱线图绘制
#SDC1
#data1<-data[which(data$gene == 'SDC1'),]
p = ggboxplot(data, x="gene", y="expression", fill = "Group",
            ylab = "Gene expression",
            xlab = "", 
            palette = c("#00A087B2","#E64B35FF"),#"#00A087B2","#E64B35FF"
            #add='jitter',
            #notch = TRUE
            )+
  theme(text = element_text(size = 10))
#p = p+rotate_x_text(60)
p=p+stat_compare_means(aes(group=Group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
p

ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 3, height = 3.5)




