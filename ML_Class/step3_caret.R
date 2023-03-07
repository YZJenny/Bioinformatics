rm(list=ls())
library(e1071)
library(mlr)
library(caret)
library(ROCR)
library(pROC)
#names(getModelInfo())

setwd('/mdshare/node9/yanzijun/Extra/230228_Model/')
data.lst <- readRDS('input_method1/data.lst')
train_EXP <- as.data.frame(t(data.lst$train.exp))
test_EXP <- as.data.frame(t(data.lst$test.exp))
EXP <- as.data.frame(rbind(train_EXP,test_EXP))
dim(EXP)

output <- read.table("Fig_method1/allDiff.txt",header = T,sep='\t') 
sigG <- output$gene[which(output$P.Val < 0.05 & abs(output$logFC) > log(1,2))]
print(length(sigG))
EXP <- EXP[,sigG]
dim(EXP)

LABEL <- rbind(data.lst$train.label,data.lst$test.label)
LABEL$label <- 'Metastasis'
LABEL$label[LABEL$Group==0] <- 'noMetastasis'
LABEL$label <- factor(LABEL$label,levels = c('Metastasis','noMetastasis'))


## 1. 数据预处理
#删除方差为0的变量
zerovar=nearZeroVar(EXP)
newdata1=EXP[,-zerovar]
dim(newdata1)

#删除强相关的变量
descrCorr = cor(newdata1)
highCorr = findCorrelation(descrCorr, 0.90)
newdata2 = newdata1[, -highCorr]
dim(newdata2)

#删除多重共线性
comboInfo = findLinearCombos(newdata2)
newdata2=newdata2[, -comboInfo$remove]
dim(newdata2)

#进行标准化并补足缺失值
Process = preProcess(newdata2)
newdata3 = predict(Process, newdata2)
dim(newdata2)

## 2.特征选择
set.seed(123)
ctrl= rfeControl(functions = rfFuncs, method = "repeatedcv",verbose = FALSE, returnResamp = "final")
Profile = rfe(newdata3, LABEL$label,rfeControl = ctrl)
print(Profile)
p1 <- plot(Profile)

## 3.数据建模及预测
newdata4=newdata3[,Profile$optVariables]
dim(newdata4)
trainx = newdata4[rownames(train_EXP),]
testx = newdata4[rownames(test_EXP),]
trainy = LABEL[rownames(train_EXP),'label']
testy = LABEL[rownames(test_EXP),'label']
#作图查看前6个变量的分布情况
p2 <- featurePlot(trainx[,1:6],trainy,plot='box')

#定义模型训练参数
set.seed(123)
fitControl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                          returnResamp = "final",classProbs = TRUE,allowParallel = T)
# fitControl <- trainControl(method = "repeatedcv", number = 10,repeats = 3,
#                            verboseIter = T,returnData = F,p = 0.75,classProbs = T, 
#                            summaryFunction = twoClassSummary,allowParallel = T)

#model1: GBM
set.seed(123)
gbmGrid = expand.grid(.n.trees = c(50, 100, 150, 200, 250, 300),.interaction.depth = c(1, 3),
                      .shrinkage = 0.1,.n.minobsinnode=10)
gbmFit = train(trainx,trainy,method = "gbm",trControl = fitControl,
                tuneGrid = gbmGrid,verbose = FALSE)
p3 <- plot(gbmFit)

#model2: random forest
set.seed(123)
rfGrid = expand.grid(mtry = c(1, 3, 5))
rfFit= train(trainx, trainy,method = "rf",trControl = fitControl,
               tuneGrid = rfGrid,verbose = FALSE)
p4 <- plot(rfFit)

#model3: svm
set.seed(123)
svmGrid = expand.grid(sigma= 2^c(-25, -20, -15,-10, -5, 0), C= 2^c(0:5))
svmFit= train(trainx, trainy,method = "svmRadial",trControl = fitControl,
             tuneGrid = svmGrid,verbose = FALSE)
p5 <- plot(svmFit)

#整合不同算法的结果
models = list(GBM = gbmFit, RF = rfFit,SVM = svmFit)
predValues = extractPrediction(models,testX = testx, testY = testy)
head(predValues)
#提取检验样本的预测结果
testValues = subset(predValues, dataType == "Test")

#计算预测概率，则使用extractProb函数
probValues = extractProb(models,testX = testx, testY = testy)
testProbs = subset(probValues, dataType == "Test")

### 4. 可视化ROC
for(m in unique(probValues$model)){
  Pred = subset(testValues, model == m)
  #观察预测结果的混淆矩阵
  confusionMatrix(Pred$pred, Pred$obs)
  prob = subset(testProbs, model == m)
  prob$lable=ifelse(prob$obs=='Metastasis',yes=1,0)
  rocobj <- roc(prob$lable, prob$Metastasis)
  print(rocobj$auc)
  # 绘制roc曲线
  
  pdf(paste('Fig_method1/',m,'_ROC_method1.pdf',sep=''),width = 5,height = 5)
  plot(rocobj,legacy.axes = TRUE,
       main=m,
       thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
       print.thres="best") # 在roc曲线上显示最佳阈值点
  text(0.8,0.9, paste("AUC = ",round(rocobj$auc,3),sep=''))
  dev.off()
}

#保存中间结果
pdf('Fig_method1/tmp_method1.pdf')
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()
save.image('Fig_method1/method1.RData')