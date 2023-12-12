library(tidyverse)
library(glmnet)
source('D://R//R-4.1.2//library//msvmRFE.R')
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest) 

setwd("E://曦昀项目//2.双疾病//1.BRCA//5.lasso筛选差异基因")
expFile="BRCA_7comm.txt"
train <- read.table(expFile, header=T, sep="\t", check.names=F)#后面svmRFE函数要求group的类型为factor
rownames(train)<-train[,1]
train<-train[,-1]
dim(train)

train[1:4,1:4]

# 转为lasso需要的格式

x <- as.matrix(train[,-1])

(y <- ifelse(train$group == "tumor", 0,1)) #把分组信息换成01

fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)

# 绘制LASSO回归曲线图

pdf("7-1A_lasso.pdf", width = 30, height = 15)

plot(fit, xvar = "dev", label = TRUE)

dev.off()

#绘制LASSO回归10折交叉验证图

cvfit = cv.glmnet(x, y,
                  
                  nfold=10, #例文描述：10-fold cross-validation
                  
                  family = "binomial", type.measure = "class")

pdf("7-2cvfit.pdf")

plot(cvfit)

dev.off()

#查看最佳lambda

cvfit$lambda.min

# 获取LASSO选出来的特征

myCoefs <- coef(cvfit, s="lambda.min")

lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]

(lasso_fea <- lasso_fea[-1])

# 把lasso找到的特征保存到文件

write.csv(lasso_fea,"7-3feature_lasso.csv") 



#SVM-REF算法输入数据
##注意一下这块的group需要是因子的形式
train$group=as.factor(train$group)
input <- train

#采用五折交叉验证 (k-fold crossValidation）

svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数

nfold = 5

nrows = nrow(input)

folds = rep(1:nfold, len=nrows)[sample(nrows)]

folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择

top.features = WriteFeatures(results, input, save=F) #查看主要变量

head(top.features)

#把SVM-REF找到的特征保存到文件

write.csv(top.features,"7-4feature_svm.csv")

# 选前300个变量进行SVM模型构建，选取的变量越多，运行速度越慢！

featsweep = lapply(1:7, FeatSweep.wrap, results, input) #300个变量

#save(featsweep,file = "featsweep.RData")

# 画图

no.info = min(prop.table(table(input[,1])))

errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#绘制基于SVM-REF算法的错误率曲线图

pdf("7-5B_svm-error.pdf",width = 5,height = 5)

PlotErrors(errors, no.info=no.info) #查看错误率

dev.off()

#绘制基于SVM-REF算法的正确率曲线图

#dev.new(width=4, height=4, bg='white')

pdf("7-6B_svm-accuracy.pdf",width = 5,height = 5)

Plotaccuracy(1-errors,no.info=no.info) #查看准确率

dev.off()

# 图中红色圆圈所在的位置，即错误率最低点

which.min(errors)

#比较lasso和SVM-REF方法一找出的特征变量，画Venn图

(myoverlap <- intersect(lasso_fea, top.features[1:which.min(errors), "FeatureName"])) #交集

#统计交叉基因有多少个

summary(lasso_fea%in%top.features[1:which.min(errors), "FeatureName"])

#绘制venn图

pdf("7C_lasso_SVM_venn.pdf", width = 15, height = 8)

grid.newpage()

venn.plot <- venn.diagram(list(LASSO = lasso_fea, #画图
                               
                               SVM_RFE = as.character(top.features[1:which.min(errors),"FeatureName"])), NULL,
                          
                          fill = c("#E31A1C","#E7B800"),
                          
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3,
                          
                          category.names = c("LASSO", "SVM_RFE"),
                          
                          main = "Overlap")

grid.draw(venn.plot)

dev.off() 