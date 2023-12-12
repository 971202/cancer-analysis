#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")
##28immune下载地址：
##http://cis.hku.hk/TISIDB/data/download/CellReports.txt
##下载完后直接txt转成gmt

##ssgsea其实也是一种免疫浸润的一种分析方法，根据表达情况和表达数据对免疫细胞进行富集分析进行

##然后将富集的情况标准到【0,1】之间

#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="GEOmatrix.txt"               #表达数据文件
#clusterFile="Subtype_Immune_Model_Based.txt"      #分型结果文件
gmtFile="28immune.gmt"              #免疫基因集文件
setwd("E://曦昀项目//2.双疾病//1.BRCA//6.ssGSEA分析")     #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,] #选取每行结果的平均数大于0的数值

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
#默认情况下，kcdf="Gaussian" 适用于输入表达式值是连续的，例如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs 或 log-TPMs。 
#当输入表达式值是整数计数时，例如来自 RNA-seq 实验的那些，那么这个参数应该设置为 kcdf="Poisson"。
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)#
#对ssGSEA打分进行矫正

##根据基因在不同样本中的表达数据的情况，得到相关的免疫细胞在不同样本中的表达的情况
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgseaScore=normalize(ssgseaScore)

#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)



##可视化
##在ssGSEA中大多数情况下只对于肿瘤样本进行计算，因此我们要先筛掉那些正常样本

#删掉正常样品
ssgseaScore=t(ssgseaScore)
group=sapply(strsplit(rownames(ssgseaScore),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
ssgseaScore=ssgseaScore[group==0,,drop=F]
rownames(ssgseaScore)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(ssgseaScore))
ssgseaScore=avereps(ssgseaScore)


##聚类信息
h_cluster = hclust(dist(t(ssgseaOut)))#计算样品距离并聚类
cluster_result=cutree(h_cluster,3)#将聚类结果分成低，中，高免疫三组
library(pheatmap)
norm_ssGSEA_Score1=read.table("ssGSEA.result.txt",sep="\t",header=T,row.names=1,check.names=F)#读取结果文件
V1<-names(cluster_result)
V2<-as.integer(cluster_result)
category<-data.frame(V1,V2)
category$V1<-as.character(category$V1)
category[,2]=paste0("Cluster",category[,2])
category=category[order(category[,2]),]
norm_ssGSEA_Score1=norm_ssGSEA_Score1[,as.vector(category[,1])]
cluster=as.data.frame(category[,2])
row.names(cluster)=category[,1]
colnames(cluster)="Cluster"
pdf("heatmap.pdf",height=5,width=9)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(norm_ssGSEA_Score1, annotation=cluster,        
         color =colorRampPalette(color.key)(50),        
         cluster_cols =F,        
         fontsize=8,        
         fontsize_row=8,        
         scale="row",        
         show_colnames=T,        
         fontsize_col=3,        
         main = "Example heatmap")
dev.off()



##根据分组进行箱线图绘制
group<-read.table("group.txt", header=T, sep="\t", check.names=F)
row.names(group) <- group[,1]
mygroup <- group[,-1]
mygroup<-as.factor(mygroup)
pdf("ssGSEA_bocplot.pdf")
draw_boxplot(ssgseaScore,mygroup,color = c("#1d4a9b","#e5171a"))
dev.off()



##ggplot画箱线图
group<-read.table("group.txt", header=T, sep="\t", check.names=F)
row.names(group) <- group[,1]
ssgseaScore=t(ssgseaScore)
rt<-cbind(group,ssgseaScore)
rt<-rt[,-1]

data=reshape2::melt(rt,id.vars=c("group"))
colnames(data)=c("group", "Immune", "Expression")

#绘制箱线图
#group=levels(factor(data$group))
data$group=factor(data$group, levels=c("normal","tumor"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="group",
                  xlab="",
                  ylab="Fraction",
                  legend.title="group",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(65)+
  stat_compare_means(aes(group=group),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

#输出图形
pdf(file="immune.diff.pdf", width=13, height=8.1)
print(boxplot)
dev.off()






