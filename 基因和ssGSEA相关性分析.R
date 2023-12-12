
setwd("E://曦昀项目//2.双疾病//1.BRCA//7.多基因和ssGSEA相关性分析")

gsea<-read.table("ssGSEA.result.txt",sep="\t",header=T,row.names=1,check.names=F)

rt<-read.table("GEOmatrix.txt",sep="\t",header=T,row.names=1,check.names=F)

gene=read.table("7commgene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(rt))
geneExp=rt[sameGene,]

#相关性检验
outTab=data.frame()


#按基因进行类型循环
for(i in rownames(geneExp)){
  exp1=geneExp[i,]
  y=as.numeric(exp1)##基因的名字的表达量
  
  #按照免疫细胞进行循环分析
  for(j in rownames(gsea)){
    outVector=data.frame(i,j)
    x=as.numeric(gsea[j,])##代表的是stromal和immune和estimate打分情况
    df1=as.data.frame(cbind(x,y))##将打分和基因的表达量放在一起
    corT=cor.test(x,y,method="spearman")##然后进行相关性分析
    cor=corT$estimate##相关性
    pValue=corT$p.value##P值
    outVector=cbind(outVector,pValue,cor)##然后就会得到肿瘤类型和基因和p值和cor值
    #p1=ggplot(df1, aes(x, y)) + 
    # xlab(j)+ylab("DNAss")+
    #ggtitle(paste0("Cancer: ",i))+theme(title=element_text(size=10))+
    #geom_point()+ geom_smooth(method="lm") + theme_bw()+
    #stat_cor(method = 'spearman', aes(x =x, y =y))
    #p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
    #if(pValue<pFilter){
    # pdf(file=paste0("DNAss.",i,"_",j,".pdf"),width=5,height=5)
    #print(p2)
    #}
    outTab=rbind(outTab,outVector)
  }
  
}
colNames=c("Gene","Immune","PValue","ssCor")
colnames(outTab)=colNames
write.table(outTab,file="ssgsea_gene.cor.txt",sep="\t",row.names=F,quote=F)


library(ggcorrplot)
library(data.table)
a<-outTab[,c("Gene","Immune","ssCor")]
genecor<-reshape2::dcast(data=a,Immune~Gene)
rownames(genecor)<-genecor$Immune
genecor<-genecor[2:ncol(genecor)]
b<-outTab[,c("Gene","Immune","PValue")]
genepvalue<-reshape2::dcast(data=b,Immune~Gene)
rownames(genepvalue)<-genepvalue$Immune
genepvalue<-genepvalue[2:ncol(genepvalue)]
pdf(file="cor.pdf", width=6, height=6)

ggcorrplot(genecor, show.legend = T, 
           p.mat = genepvalue, digits = 2, tl.cex=8, 
           sig.level = 0.05,insig = 'blank',lab = F)
dev.off()

outTab$pstar <- ifelse(outTab$PValue < 0.05,
                       ifelse(outTab$PValue  < 0.01,
                              ifelse(outTab$PValue  < 0.001,"***","**"),"*"),"")


pdf(file="cortest.pdf", width=10, height=6)                                    
ggplot(outTab, aes(Immune, Gene)) +
  geom_tile(aes(fill = ssCor), colour = "white",size=1) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  geom_text(aes(label=pstar),col ="black",size = 5) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0("* p < 0.05","\n\n",
                    "** p < 0.01","\n\n",
                    "*** p < 0.001","\n\n",
                    "Correlation"))

dev.off()