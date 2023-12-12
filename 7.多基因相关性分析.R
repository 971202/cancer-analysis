
setwd("E://曦昀项目//2.双疾病//1.BRCA//8.基因相关性")

rt<-read.table("GEOmatrix.txt",sep="\t",header=T,row.names=1,check.names=F)

gene=read.table("7commgene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(rt))
geneExp=rt[sameGene,]
a<-t(geneExp)
##相关性分析
corr<-cor(a)
library(corrplot)
pdf(file = "cor.pdf")
corrplot(corr =corr,method = "circle",type = "upper", tl.pos="lt",insig="blank",sig.level = 0.1, pch.cex = 0.7) # 默认p>0.05的相关性系数不展示。设置sig.level = 0.1表示表示p>0.1的相关性系数都不展示。
corrplot(corr = corr,type="lower",add=TRUE,method="number", tl.pos = "n",cl.pos = "n",diag=FALSE,col = "black" ) 

dev.off()

