
setwd("E://曦昀项目//2.双疾病//1.BRCA//9.roc诊断差异基因")
rt<-read.table("GSE28242.txt",sep="\t",header=T,row.names=1,check.names=F)

gene=read.table("7commgene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(rt))
geneExp=rt[sameGene,]
a<-t(geneExp)
group<-read.table("group.txt", header=T, sep="\t", check.names=F)
rownames(group)<-group$ID
brca<-cbind(group,a)
brca<-brca[,-1]
brca<-brca[order(brca$group),]
brca$group<-ifelse(brca$group=="normal","0","1")
brca$group<-as.numeric(brca$group)

roc_FCER1G<-roc(group~FCER1G,data=brca,levels=c("0","1"))  
roc_FCER1G
roc_ITGAM<-roc(group~ITGAM,data=brca,levels=c("0","1")) 
roc_ITGAM
roc_LCP2<-roc(group~LCP2,data=brca,levels=c("0","1"))  
roc_LCP2
roc_LILRB2<-roc(group~LILRB2,data=brca,levels=c("0","1"))  
roc_LILRB2
roc_MNDA<-roc(group~MNDA,data=brca,levels=c("0","1"))  
roc_MNDA
roc_SPI1<-roc(group~SPI1,data=brca,levels=c("0","1"))  
roc_SPI1
roc_TYROBP<-roc(group~TYROBP,data=brca,levels=c("0","1"))  
roc_TYROBP


##绘制曲线
plot(smooth(roc_FCER1G),print.auc=TRUE,           #去掉smooth为不光滑曲线
                   print.auc.x=0.4,print.auc.y=0.7,
                   print.auc.col=c("indianred2"),
                   auc.polygon=TRUE,
                   auc.polygon.col="linen",  #阴影颜色 #E16A86  #00AD9A  #FFFFCC  beige seashell
                    grid=c(0.5,0.2),
                   grid.col=c("black","black"),
                   print.thres=FALSE, 
                   main="ROC",
               col="indianred2",
                   legacy.axes=TRUE)    #拖曳窗口改变坐标轴
#再添加一条曲线（绘制多条曲线于同一面板）
plot(roc_ITGAM,add=TRUE,col = "dodgerblue3",print.thres = FALSE,
               print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.6)  #颜色  00AD9A  steelblue dodgerblue2
plot(roc_LCP2,add=TRUE,col = "#00AD9A",print.thres = FALSE,
           print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.5)
plot(roc_LILRB2,add=TRUE,col = "#E16A86",print.thres = FALSE,
     print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.4)
plot(roc_MNDA,add=TRUE,col = "#ff4500",print.thres = FALSE,
     print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.3)
plot(roc_SPI1,add=TRUE,col = "#006666",print.thres = FALSE,
     print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.2)
plot(roc_TYROBP,add=TRUE,col = "#00cc33",print.thres = FALSE,
     print.auc=TRUE,print.auc.x=0.4,print.auc.y=0.1)







##多变量分析
glmfit <- glm(group ~ C2CD4A+S100A7+TIGD1 , data =  brca)  #group是数值型
print(glmfit)
predict(glmfit)
brca$predicted_values <- predict(glmfit)  #将结果加到表格中
head(brca)
fitted(glmfit)    #结果与predict(glmfit)或predict(glmfit,type = "response")相同
exp(cbind(OR=coef(glmfit),confint(glmfit))) #查看OR及置信区间
roc_rbp_penk_gcgr_inha <- roc(brca$group,fitted(glmfit),plot=TRUE) #AMI$status是因变量，fitted(glmfit)是自变量
pdf(file="ROCtogh.pdf", width=5, height=5)
plot(roc_rbp_penk_gcgr_inha,print.auc=TRUE,           #去掉smooth为不光滑曲线
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.col=c("tan4"),
     auc.polygon=TRUE,
     auc.polygon.col="beige",  #阴影颜色 #E16A86  #00AD9A  #FFFFCC  beige seashell  linen
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     print.thres=TRUE, 
     main="ROC",
     col="tan4",
     legacy.axes=TRUE) 

dev.off()

library(pROC)
res<-roc(group~.,data=brca,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=F,# 是否平滑曲线
         levels=c('tumor','normal'),direction=">" #设置分组方向
)

library(ggplot2)
pdf("new_roc.pdf",width = 10,height =8)
p<- ggroc(res, legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw()+  # 设置背景
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid = element_blank())

p+annotate("text",x=0.85,y=0.055,label=paste("FCER1G-AUC = ", round(res$FCER1G$auc,3)))+
  annotate("text",x=0.85,y=0.15,label=paste("ITGAM-AUC = ", round(res$ITGAM$auc,3)))+
  annotate("text",x=0.85,y=0.25,label=paste("LCP2-AUC = ", round(res$LCP2$auc,3)))+
  annotate("text",x=0.85,y=0.35,label=paste("LILRB2-AUC = ", round(res$LILRB2$auc,3)))+
  annotate("text",x=0.85,y=0.45,label=paste("MNDA-AUC = ", round(res$MNDA$auc,3)))+
  annotate("text",x=0.85,y=0.55,label=paste("SPI1-AUC = ", round(res$SPI1$auc,3)))+
  annotate("text",x=0.85,y=0.65,label=paste("TYROBP-AUC = ", round(res$TYROBP$auc,3)))
#plegend(col=c("#FA8072","#FFD700","#32CD32"), lwd=2, bty = 'n')
dev.off()
