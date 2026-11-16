

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\03genesets\\ARGs'
setwd(path)

#差异基因与切入点基因集取交集
genes<-read.table('ARGs.txt',sep='\t')
diff<-read.table('DEGs.txt',sep='\t')
wgcna<-read.table('wgcna_hubs.txt',sep='\t')

a<-genes$V1
b<-diff$V1
length(a);length(b)
#[1] 1542
#[1] 779

d<-wgcna$V1
length(d)
#[1] 3267

c<-intersect(a,b)
length(c)
#[1] 93
c<-intersect(c,d)
length(c)
#[1] 1

#save
write.table(data.frame(c),'diff_intergene.txt',sep='\t',
     col.names=F,row.names=F,quote=F)
#write.table(data.frame(c),'wgcna_intergene.txt',sep='\t',
#     col.names=F,row.names=F,quote=F)


#####交集分析绘图1
#upset
library(UpSetR)

p<-list('ARGs'=a,'DEGs'=b)

#save plot
pdf("upset.pdf", onefile=F, width =5, height =4)
upset(fromList(p), sets.bar.color=c("red", "blue"),
     text.scale = c(2, 1.5, 1.5, 1.5, 2, 2))
dev.off()

#####交集分析绘图2
library(VennDiagram)

s1 <- list('HMGs'=a,'DEGs'=b)

v1 <- venn.diagram(x = s1, 
                   filename = NULL, # 直接给一个名称会自动保存文件到本地
                   # 输出的图形参数
                   # imagetype = "png", # 输出图片类型(tiff,png,svg)
                   # height = 1000, # 图片高度
                   # width = 1000, # 图片宽度
                   # resolution = 300, # 图片分辨率
                   
                   scaled = F, # 根据比例显示大小
                   alpha=c(0.3, 0.3), # 设置每个区块的透明度
                   
                   ## 下面是除了标题外，图形其他元素的设置参数
                   # 图形元素设置：圈
                   lwd = 1, # 圆圈线条的粗细：1 2 3 4 5 6
                   lty = 'blank', # 圆圈线条的类型：1为实线，2为虚线，blank为无线条
                   col = c("black","red"), # 圆圈线条颜色
                   fill = c("skyblue", "springgreen"),  # 圆圈颜色
                   
                   # 图形元素设置：数字
                   cex = 1, # 数字大小
                   fontface = "bold", # 加粗
                   fonrfamily = "sans", # 数字字体
                   
                   
                   # 图形元素设置：标签即(category)
                   cat.cex = 1,  # 标签字体大小
                   cat.col = "black",  # 标签字体色
                   cat.fontface = "bold",  # 加粗
                   cat.default.pos = "outer",  # 标签内外位置, 在圆圈内还是圆圈外，outer 内 text 外
                   cat.pos = c(0, 0),  # 标签旋转位置，用圆的度数
                   cat.dist = c(0.05,0.05),  # 标签离圆圈位置,离圆的距离，如果标签与圆圈重叠，可以调整这个参数
                   cat.fontfamily = "sans",  # 标签字体
)

cowplot::plot_grid(v1)

pdf("venn.pdf", onefile=F, width =4, height =4)
cowplot::plot_grid(v1)
dev.off()





