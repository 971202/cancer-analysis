##fpkm值进行差异分析

##为什么不用fpkm值进行差异分析，因为fpkm值考虑到了基因的长度对基因的表达量的影响，所以才进行了归一化的影响（同一个样本比较基因的表达量

##但是在进行差异分析的时候在不同的样本中基因的表达长度不同，也不考虑这个因素，所以用这个没有意思

##fpkm做差异分析的时候，可以比较同一个样本中所有的基因的差异

##比如热图信息，就是比较所有的基因在同一个样本中的表达情况如果，所以在画热图的时候，我们使用了log(data+1)进行

##RPKM和FPKM值的区别，RPKM是针对于单端进行，FPKM是针对于双端进行，所以这两个其实是一个道理
##TPM其实是RPKM的百分比
##标准化程度而言，FPKM值是不如TPM值的

##安装软件包
# install.packages("pacman")          #安装软件包
library(pacman)                      
p_load(limma,edgeR,pheatmap)   

##读取数据
FCfileter=1
pfilter=0.05
setwd("E://曦昀项目//2.双疾病//1.BRCA//2.差异分析")

expfile="GEOmatrix.txt"
a=read.table(expfile,sep='\t',quote = "",fill = T,
             comment.char = "!",header = T) # 提取表达矩阵
rownames(a)=a[,1]
a <- a[,-1]
##将fpkm值转化成tpm值
expMatrix <- a
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
##这个信息需要对数据进行更改，因为根据不同的数据有不同的组别
group_list=c(rep("tumor",10),rep("normal",6))
## 强制限定顺序
group_list <- factor(group_list,levels = c("tumor","normal"),ordered = F)


#表达矩阵数据校正
exprSet <- tpms

pdf(file="before_normalization.pdf", width=10, height=7)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2,xaxt="s",yaxt="n",cex.axis=0.8,main="Before Normalization")
dev.off()
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
pdf(file="normalization.pdf", width=10, height=7)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2,xaxt="s",yaxt="n",cex.axis=0.8,main="Normalization")
dev.off()
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
group_list=c(rep("tumor",10),rep("normal",6))

design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')

bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
#bp(dat['BAT2L1',])
#ggsave("BAT2L1_box.png")
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
#save(deg,file = 'deg.Rdata')
write.table(deg,"All_gene.txt",sep="\t",quote = F)
Diffgene = deg[(deg$adj.P.Val < pfilter & (deg$logFC>=FCfileter | deg$logFC<=(-FCfileter))),]
write.table(Diffgene, "Diff_gene05.txt",sep="\t",quote=F)
Upgene = deg[(deg$adj.P.Val < pfilter & (deg$logFC>=FCfileter)),]
write.table(Upgene, "Up_gene05.txt",sep="\t",quote=F)
Downgene = deg[(deg$adj.P.Val < pfilter & (deg$logFC<=(-FCfileter))),]
write.table(Downgene, "Down_gene05.txt",sep="\t",quote=F)


#热图的构建
heatExp<-dat[rownames(Diffgene),]
#heatExp=log2(heatExp+1)
heatExp=heatExp[1:30,]
Type=c(rep("T",10),rep("N",6))   
names(Type)=colnames(heatExp)
Type=as.data.frame(Type)
ann_colors = list(
  group = c(low="#32CD32", high="#DC143C"))

pdf("heatmap.pdf",20,10)
pheatmap(heatExp, 
         #scale="row", #妯悜鏍囧噯鍖栨瘮杈?
         annotation_col=Type, 
         color = colorRampPalette(c("#1E90FF", "white", "#DC143C"))(50),
         cluster_cols =F,
         fontsize = 10,fontsize_row=10,fontsize_col=5,
         show_colnames = F,
         annotation_legend = T,
         annotation_names_col = T,
         annotation_colors =ann_colors[1]
)
dev.off()

#缁樺埗鐏北鍥?
canoExp<-deg[abs(deg$logFC)< FCfileter | deg$adj.P.Val>=pfilter,]
xmax<-max(deg$logFC)
ymax<-max(-log10(deg$adj.P.Val))
downgene<- transform(Downgene,adj.P.Val=-log10(Downgene$adj.P.Val))
upgene<- transform(Upgene,adj.P.Val=-log10(Upgene$adj.P.Val))
canoExp<- transform(canoExp,adj.P.Val=-log10(canoExp$adj.P.Val))
pdf("Volcano.pdf")
plot(canoExp$logFC,canoExp$adj.P.Val,xlim = c(-xmax,xmax),ylim=c(0,ymax),col="black",
     pch=16,cex=0.9,main = "Volcano",xlab = "logFC",ylab="-log10(adj.P.Val)")
points(upgene$logFC,upgene$adj.P.Val,col="#DC143C",pch=16,cex=0.9)
points(downgene$logFC,downgene$adj.P.Val,col="#1E90FF",pch=16,cex=0.9)
abline(v=0,lwd=3,lty=2)
dev.off()




##和上面的差异基因分析的第二种方法，下面的可以不运行，两个代码生成一样的结果


library(clusterProfiler)
library(enrichplot)
library(tidyverse)

## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
if(T){
  logFC_t=0.585
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  
  save(DEG,file = 'anno_DEG.Rdata')
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
}

##火山图
## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable',
              ifelse( df$logFC >1.5,'up',
                      ifelse( df$logFC < -1.5,'down','stable') )
  )
  table(df$g)#此时获得的上调基因为20，下调基因为145
  df$symbol=rownames(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "symbol", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = rownames(head(head(deg))),
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  
  
}
ggsave("volcano.png")


##热图的绘制

if(T){ 
  #load(file = 'step2-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  #load(file = 'deg.Rdata')
  x=deg$logFC
  names(x)=rownames(deg)
  cg=c(names(head(sort(x),20)),#老大的代码中是100、100，由于我获得上下调基因少，就改了 
       names(tail(sort(x),145)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F)
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F,
           annotation_col=ac)
  
  
}

#ggsave("heatmap.png")

