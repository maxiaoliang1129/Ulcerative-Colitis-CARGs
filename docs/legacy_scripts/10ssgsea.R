

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\10ssgsea'
setwd(path)

#data
tpm<-read.table('normalize.txt',header=T,row.names=1,sep='\t')

expr_tpm<-tpm
dim(expr_tpm)
#[1] 39374   1320


#导入所需包
library(org.Hs.eg.db)
library(stringr)
library(BiocGenerics)
library(clusterProfiler)
library(enrichplot)
library(future)
library(future.apply)

#计算相关性
#定义函数:获得目标基因与其他基因之间的P值和Cor值
exprSet<-expr_tpm

batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="pearson")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

#读入特征基因
gene<-read.table('final.model.txt',header = F,sep='\t')
gene
#       V1
#1 CCL11
#2 CLDN1


######################单基因gsea1:CLDN1#################
dd<-batch_cor("CLDN1")
dim(dd)

#CCL11的相关性阈值是0.05和0.5 CLDN1的相关性阈值设置为0.85
dd<-dd[which(dd$p.value < 0.05),]
dd<-dd[which(abs(dd$cor) > 0.5),]
dim(dd)
#[1] 137   4  (0.5)
#[1] 764    4  (0.5)

gene <- dd$mRNAs

#换成大写
gene<-toupper(gene)
dd$mRNAs<-toupper(dd$mRNAs)

#转化格式
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_df <- data.frame(cor=dd$cor,SYMBOL = dd$mRNAs)

#去重复
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
geneList <- gene_df$cor
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T)

dotplot_internal <- function(object, x = "GeneRatio", color = "pvalue",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x", decreasing=TRUE) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) + scale_size(range=c(3, 8))
  
}


##'@GSEA GO 富集
##'读取go.gmt文件
#GOgmt<-read.gmt("h.all.v2024.1.Hs.symbols.gmt")
GOgmt<-read.gmt("c5.go.bp.v2024.1.Hs.symbols.gmt")
##'@富集分析
GO<-GSEA(geneList,TERM2GENE = GOgmt)
go_result<-GO@result

#save
write.table(go_result,'CLDN1_go.txt',col.names = T,row.names = T,sep='\t',quote=F)

#plot
#先排序
go_result<-go_result[order(go_result$pvalue),]
gseaplot2(GO,1:3,color="red",pvalue_table = FALSE)

#save
pdf(file="CLDN1_go.pdf", width = 5, height = 4)
gseaplot2(GO,1:3,color="red",pvalue_table = FALSE)
dev.off()

