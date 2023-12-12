#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")
##28immune���ص�ַ��
##http://cis.hku.hk/TISIDB/data/download/CellReports.txt
##�������ֱ��txtת��gmt

##ssgsea��ʵҲ��һ�����߽����һ�ַ������������ݱ�������ͱ������ݶ�����ϸ�����и�����������

##Ȼ�󽫸����������׼����0,1��֮��

#���ð�
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="GEOmatrix.txt"               #���������ļ�
#clusterFile="Subtype_Immune_Model_Based.txt"      #���ͽ���ļ�
gmtFile="28immune.gmt"              #���߻����ļ�
setwd("E://������Ŀ//2.˫����//1.BRCA//6.ssGSEA����")     #���ù���Ŀ¼

#��ȡ���������ļ�,���������ļ�����
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,] #ѡȡÿ�н����ƽ��������0����ֵ

#��ȡ�����ļ�
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA����
#Ĭ������£�kcdf="Gaussian" �������������ʽֵ�������ģ���������߶ȵ�΢����ӫ�ⵥԪ��RNA-seq log-CPMs��log-RPKMs �� log-TPMs�� 
#���������ʽֵ����������ʱ���������� RNA-seq ʵ�����Щ����ô�������Ӧ������Ϊ kcdf="Poisson"��
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)#
#��ssGSEA��ֽ��н���

##���ݻ����ڲ�ͬ�����еı������ݵ�������õ���ص�����ϸ���ڲ�ͬ�����еı�������
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#�����ǽ�ÿ�������в�ͬ������ϸ��������׼����0-1֮��
ssgseaScore=normalize(ssgseaScore)

#���ssGSEA��ֽ��
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)



##���ӻ�
##��ssGSEA�д���������ֻ���������������м��㣬�������Ҫ��ɸ����Щ��������

#ɾ��������Ʒ
ssgseaScore=t(ssgseaScore)
group=sapply(strsplit(rownames(ssgseaScore),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
ssgseaScore=ssgseaScore[group==0,,drop=F]
rownames(ssgseaScore)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(ssgseaScore))
ssgseaScore=avereps(ssgseaScore)


##������Ϣ
h_cluster = hclust(dist(t(ssgseaOut)))#������Ʒ���벢����
cluster_result=cutree(h_cluster,3)#���������ֳɵͣ��У�����������
library(pheatmap)
norm_ssGSEA_Score1=read.table("ssGSEA.result.txt",sep="\t",header=T,row.names=1,check.names=F)#��ȡ����ļ�
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



##���ݷ����������ͼ����
group<-read.table("group.txt", header=T, sep="\t", check.names=F)
row.names(group) <- group[,1]
mygroup <- group[,-1]
mygroup<-as.factor(mygroup)
pdf("ssGSEA_bocplot.pdf")
draw_boxplot(ssgseaScore,mygroup,color = c("#1d4a9b","#e5171a"))
dev.off()



##ggplot������ͼ
group<-read.table("group.txt", header=T, sep="\t", check.names=F)
row.names(group) <- group[,1]
ssgseaScore=t(ssgseaScore)
rt<-cbind(group,ssgseaScore)
rt<-rt[,-1]

data=reshape2::melt(rt,id.vars=c("group"))
colnames(data)=c("group", "Immune", "Expression")

#��������ͼ
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

#���ͼ��
pdf(file="immune.diff.pdf", width=13, height=8.1)
print(boxplot)
dev.off()





