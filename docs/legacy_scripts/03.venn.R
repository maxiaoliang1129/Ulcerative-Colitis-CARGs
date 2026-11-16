

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\03venn\\00clusterprofiler'
setwd(path)

#通路分析
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)

#查看支持转换的类型
keytypes(org.Hs.eg.db)

#全部差异基因做富集分析
gene<-read.table('diff_intergene.txt',header=F,sep='\t')
symbol<-gene$V1

length(symbol)
#[1] 93

entrez<-bitr(symbol,fromType = "SYMBOL",
             toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
dim(entrez)
#[1] 93  2

test<-entrez$ENTREZID
head(test)

##内置示例数据（确定背景基因集）
#data(geneList, package="DOSE") #富集分析的背景基因集

#这里选择所有已知的基因作为背景
background_genes <- keys(org.Hs.eg.db)

ego_ALL <- enrichGO(gene = test, 
                    universe = background_genes, #背景基因集
                    #universe = names(geneList), #背景基因集
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    minGSSize=1,
                    maxGSSize=500,
                    #pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    #qvalueCutoff = 1,
                    pvalueCutoff = 0.05, #P值会过滤掉很多，0.05
                    qvalueCutoff = 0.1, #fdr值 0.05 0.01 0.1（样本较少）
                    readable = TRUE) #Gene ID 转成gene Symbol ，易读

GO_result<-ego_ALL@result
dim(GO_result)
#[1] 506  13

write.csv(GO_result,"GO_p0.05result.csv",row.names =FALSE)


#KEGG分析
#plot1 p=0.05
kk<-enrichKEGG(gene = test,organism = 'hsa',
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               universe = background_genes)
dim(kk)
#[1] 37 14

write.csv(as.data.frame(kk),"KEGG_p0.05result.csv",row.names =FALSE)




rm(list=ls())

path<-'D:\\LZ-25234\\UC\\03venn\\00clusterprofiler'
setwd(path)

library(ggplot2)
library(openxlsx)

#data
#GO数据
go<-read.csv('GO_p0.05result.csv',header = T,sep = ',')
head(go);dim(go) #[1] 506  13
go_enrich<-go
#数据处理#
#go_enrich$term<-paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term<-go_enrich$Description
go_enrich$term<-factor(go_enrich$term, levels = go_enrich$term,ordered = T) #转成因子，防止重新排列

#各取top10个基因
BP<-go_enrich[which(go_enrich$ONTOLOGY == 'BP'),]
CC<-go_enrich[which(go_enrich$ONTOLOGY == 'CC'),]
MF<-go_enrich[which(go_enrich$ONTOLOGY == 'MF'),]
dim(BP);dim(CC);dim(MF)
#[1] 407  14
#[1] 43 14
#[1] 56 14

BP<-BP[order(BP$pvalue),]
CC<-CC[order(CC$pvalue),]
MF<-MF[order(MF$pvalue),]

go_enrich<-rbind(BP[1:5,],CC[1:5,],MF[1:5,])
dim(go_enrich) #[1] 15 14


###########################多种绘图方法#########################
##纵向柱状图-根据ONTOLOGY类型绘制#
p1 <- ggplot(go_enrich,
             aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')


##横向柱状图-根据ONTOLOGY类型绘制#
p2 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  #labs(title = "GO enrichment analysis results")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p2 
#根据ONTOLOGY分类信息添加分组框#
p2 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

##纵向柱状图——-根据pvalue值绘制#
p3 <- ggplot(go_enrich,aes(y=term,x=Count,fill=pvalue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "",  #设置标题、x轴和Y轴名称
       x = "Gene number", 
       y = "GO Terms")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
p3
#根据ONTOLOGY分类信息添加分组框#
p3+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

##横向柱状图-根据pvalue值绘制#
p4 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=pvalue)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_gradient(low = "red",high ="blue" ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p4 
#根据GO富集分析分类信息添加分组框#
p4 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

##气泡图#
P5 <- ggplot(go_enrich,
             aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
P5
#根据ONTOLOGY分类信息添加分组框#
P5 + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

###############################################################
#保存图片
#tiff(file='GO.barplot.tiff',width = 6,height = 5,units = "in",bg = "white",res = 350)
pdf(file = "GO.barplot.pdf",height = 5,width = 7)
p1
dev.off()





rm(list=ls())

path<-'D:\\LZ-25234\\UC\\03venn\\00clusterprofiler'
setwd(path)

library(ggplot2)
library(openxlsx)

#data
#GO数据
go<-read.csv('KEGG_p0.05result.csv',header = T,sep = ',')
head(go)
go_enrich<-go
#数据处理#
#go_enrich$term<-paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term<-go_enrich$Description
go_enrich$term<-factor(go_enrich$term, levels = go_enrich$term,ordered = T) #转成因子，防止重新排列
dim(go_enrich)
#[1] 37 15

#排序并取前15
go_enrich<-go_enrich[order(go_enrich$p.adjust),]
go_enrich<-go_enrich[1:15,]

###########################多种绘图方法#########################
##纵向柱状图-根据ONTOLOGY类型绘制#
p1 <- ggplot(go_enrich,
             aes(x=term,y=Count, fill=subcategory)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("KEGG term") +  #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "")+  #设置标题
  theme_bw()
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(subcategory~., scale = 'free_y', space = 'free_y')


##横向柱状图-根据ONTOLOGY类型绘制#
p2 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=subcategory)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p2 
#根据ONTOLOGY分类信息添加分组框#
p2 + facet_grid(.~subcategory, scale = 'free_x', space = 'free_x')

##纵向柱状图——-根据pvalue值绘制#
p3 <- ggplot(go_enrich,aes(y=term,x=Count,fill=pvalue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "",  #设置标题、x轴和Y轴名称
       x = "Gene number", 
       y = "KEGG Terms")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
p3
#根据ONTOLOGY分类信息添加分组框#
p3+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

##横向柱状图-根据pvalue值绘制#
p4 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=pvalue)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_gradient(low = "red",high ="blue" ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p4 
#根据GO富集分析分类信息添加分组框#
p4 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

##气泡图#
P5 <- ggplot(go_enrich,
             aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
P5
#根据ONTOLOGY分类信息添加分组框#
P5 + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

###############################################################
#保存图片
#tiff(file='KEGG.barplot.tiff',width = 8.5,height = 5,units = "in",bg = "white",res = 350)
pdf(file = "KEGG.barplot.pdf",height = 5,width = 8.5)
p1
dev.off()




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\03venn\\01pheatmap'
setwd(path)

#热图
headmap_expr<-read.table('normalize.txt',header=T,sep='\t',row.names = 1)
expr<-headmap_expr
dim(expr)
#[1] 39374  1320

#分组信息
group<-read.table('group.txt',sep='\t',header=T)
colnames(group)<-c('Sample','Group')

#group 排序
group<-group[order(group$Group),]

rownames(group)<-group[,1]
sample<-rownames(group)
group<-data.frame(group$Group)
rownames(group)<-sample
colnames(group)[1]<-'Group'
group$Group<-as.factor(group$Group)

expr<-expr[,rownames(group)]

#只挑选细胞黏附分子相关基因
gene<-read.table('diff_intergene.txt',header = F,sep = '\t')
gene<-gene$V1

rt<-expr[gene,]
dim(rt)
#[1] 93  1320
rt<-na.omit(rt)
dim(rt) #[1] 93  1320

#修改行名
df<-rt
old_names <- rownames(df)
old_names[old_names == "PVR"] <- "Pvr"
rownames(df) <- old_names
rt<-df

#plot heatmap
library(pheatmap)

ann_color <- list(Group = c(control = "#4DBBD5FF", UC = "#DC0000FF"))
pheatmap(rt, annotation=group, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows =T,
         fontsize=8,
         fontsize_row=6,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize_col=3,
         annotation_colors = ann_color)

#save plot
pdf(file="heatmap.pdf", width = 5, height = 10,onefile = F)
pheatmap(rt, annotation=group, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows =T,
         fontsize=8,
         fontsize_row=6,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()

