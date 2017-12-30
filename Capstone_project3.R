###Load Dataset
T2D_df.T1<-read.csv("T2D_df.T3.csv",header=T)###*** Run This
### ***Data Pre-processing steps for Non-diabetic dataset1***
###load packages
require(pROC)
require(caret)
require(DMwR)
require(dplyr)
require(class)
require(R.utils)
require(glmnet)
library(glmnet)
require(magrittr)
library(magrittr)
require(ggplot2)
###set the seed to 1234
set.seed(1234)
###Data Exploration & Inspection
str(T2D_df.T1)
dim(T2D_df.T1)
names(T2D_df.T1)
class(T2D_df.T1)
T2D_df.T1$Cell.Type
is.factor(T2D_df.T1$Cell.Type)
levels(T2D_df.T1$Cell.Type)
target<-T2D_df.T1$Cell.Type
###***Check the balance of the class***
print(table(T2D_df.T1$Cell.Type))
print(prop.table(table(T2D_df.T1$Cell.Type))) ###There is class imbalance***
### Filter Data: Remove genes with low expression
T2D_df.T1$Cell.Type<-NULL

myfunction <- function(x){
  
  {    P<-mean(x)
  }
  return(P)
}
alist=list() 
blist=list()
flist=list()
for(i in 1:length(T2D_df.T1)){
  P<-myfunction(T2D_df.T1[,i])
  
  if (P >2){
    print(paste("The column is", i, "P >2" ))
    alist = list(alist, i) # combining them in a list
    
  }
  
  else  if(P<=2)
    print(paste("The column is", i, "skip and omit" ))
  blist = list(blist, i) # combining them in a list
}

desiredList<-do.call(c, list(alist))
undesiredList<-do.call(c, list(blist))

desiredVector<-unlist(desiredList, use.names=FALSE)
undesiredVector<-unlist(undesiredList, use.names=FALSE)

head(desiredVector)
tail(desiredVector)
head(undesiredVector)
tail(undesiredVector)
length(desiredVector)
length(undesiredVector)

T2D_df.T1Subset<- subset(T2D_df.T1, select = desiredVector)
dim(T2D_df.T1Subset)   ####9273
head(T2D_df.T1Subset[1:12])
tail(T2D_df.T1Subset[1:12])
head(T2D_df.T1Subset[1:4])

dfNew<-cbind(target,T2D_df.T1Subset)
dim(dfNew)
head(dfNew[1:12])
tail(dfNew[1:12])
levels(dfNew$target)
print(table(dfNew$target))
write.csv(dfNew, file = "dfNew2.csv", row.names = FALSE)
dfNew<-read.csv("dfNew2.csv",header=T)
dim(dfNew)
sink("AnovaTest5.doc")
for(j in 2:length(dfNew)){
  #### create side by side boxplot of the data:
  #boxplot(dfNew[,j]~dfNew$target)
  ###Create an ANOVA table:
  
  columns <- names(dfNew[j])
  anovaresult<- anova(aov(dfNew[,j]~target,data=dfNew)) 
  ###Conduct a Tukey's multiple comparison precedure:
  #posthocresult <-TukeyHSD(aov(dfNew[,j]~target,data=dfNew))
  if(anovaresult$Pr[1] < 0.05){ 
    flist = list(flist, j)}
  else 
    print(columns)
    

}    
sink()

length(unlist(flist))
desiredVector2<-unlist(flist, use.names=FALSE)
T2D_df.T1Subset2<- subset(dfNew, select = desiredVector2)
T2D_df_Anova<-cbind(target,T2D_df.T1Subset2)
dim(T2D_df_Anova)
head(T2D_df_Anova[1:12])
tail(T2D_df_Anova[1:12])
levels(T2D_df_Anova$target)
print(table(T2D_df_Anova$target))
###Removing replicated variables:
set.seed(1234)
dim(T2D_df_Anova)
df2 = cor(T2D_df_Anova[,-1])
hc = findCorrelation(df2, cutoff=0.99) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = T2D_df_Anova[,-1][,-c(hc)]
dim(reduced_Data)
head(reduced_Data[,1:12])
colnames(reduced_Data)
T2D_reduced_Data<-cbind(target,reduced_Data)
dim(T2D_reduced_Data) ###***5175
head(T2D_reduced_Data[,1:12])
write.csv(T2D_reduced_Data, file = "T2D_reduced_Data.csv", row.names = FALSE)
T2D_reduced_Data<-read.csv("T2D_reduced_Data.csv",header=T)

#Lasso Regression Step
set.seed(1234)
cv.glmnet.fit <- cv.glmnet(as.matrix(T2D_reduced_Data[,-1]),T2D_reduced_Data[,1],alpha = 0.5, type.multinomial ="grouped",type.measure = "class",nlambda =200, family ="multinomial" )
plot(cv.glmnet.fit)
tmp<-coef(cv.glmnet.fit,s="lambda.min")
T<-as.matrix(tmp[[1]])
coefs<-data.frame(T[(abs(T[,1]) >0),])
features<-row.names(coefs)
features
idx<-match(features,colnames(T2D_reduced_Data),nomatch = 0)
T2D_reduced_Data2<-T2D_reduced_Data[,idx]
T2D_reduced_Data2$target<-T2D_reduced_Data$target
head(T2D_reduced_Data2[90:96])
tail(T2D_reduced_Data2[90:96])
dim(T2D_reduced_Data2)
colnames(T2D_reduced_Data2)[96] <- "Class.Label"
dput(names(T2D_reduced_Data2))
colnames(T2D_reduced_Data2)<-c("V92", "V170", "V287", "V386", "V440", "V506", 
                                  "V529", "V618", "V816", "V942", "V979", "V1052", "V1389", "V1497", 
                                  "V1666", "V1785", "V1788", "V1921", "V1923", "V1969", "V2001", 
                                  "V2167", "V2252", "V2309", "V2465", "V2503", "V2514", "V2673", 
                                  "V2717", "V2928", "V3057", "V3117", "V3253", "V3260", "V3296", 
                                  "V3545", "V3712", "V3784", "V3884", "V4015", "V4050", "V4377", 
                                  "V4416", "V4553", "V4703", "V4724", "V4932", "V4965", "V5015", 
                                  "V5078", "V5165", "V5292", "V5296", "V5317", "V5448", "V5462", 
                                  "V5773", "V5804", "V6137", "V6958", "V8205", "V8223", "V8634", 
                                  "V8747", "V9100", "V9344", "V9449", "V10489", "V10563", "V11634", 
                                  "V11765", "V12116", "V12163", "V12351", "V12636", "V12658", "V12961", 
                                  "V12998", "V13063", "V13187", "V15145", "V15856", "V16869", "V16870", 
                                  "V17259", "V17914", "V18602", "V19238", "V22129", "V26473", "V26965", 
                                  "V28066", "V29326", "V34807", "V35693", "Class.Label")
write.csv(T2D_reduced_Data2, file = "T2D_reduced_Data2.csv", row.names = FALSE)
T2D_reduced_Data2<-read.csv("T2D_reduced_Data2.csv",header=T)

ddgg<-T2D_reduced_Data2 %>% count(Label = factor(Class.Label), Class.Label = factor(Class.Label)) %>% 
  ungroup() %>%
  mutate(percentage = round(prop.table(n) * 100),digits=2) %>% 
 ggplot(aes(x = Class.Label, y = percentage, fill = Class.Label)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  ylim(0, 100)+
  geom_text(aes(y = percentage + .5,    # nudge above top of bar
                label = paste0(percentage, '%'),size=14),    # prettify
            position = position_dodge(width = .9), 
            size = 5)+theme_light() +ggtitle("Imbalanced Class 'Endocrine Pancreas: Islets of Langerhans' Dataset3")
ddgg+theme(axis.text=element_text(size=12),
             axis.title=element_text(size=14,face="bold"))
###*** Training and Testing***###
dim(T2D_reduced_Data2)
T2D_reduced_Data2S<-scale(T2D_reduced_Data2[,-96])
T2D_reduced_Data2S<-as.data.frame(T2D_reduced_Data2S)
T2D_reduced_Data2S<-cbind(T2D_reduced_Data2S,T2D_reduced_Data2$Class.Label)
colnames(T2D_reduced_Data2S)[96] <- "Class.Label"
colnames(T2D_reduced_Data2S)
class(T2D_reduced_Data2S)

train<-sample_frac(T2D_reduced_Data2S, 0.70)
sid<-as.numeric(rownames(train)) # because rownames() returns character
test<-T2D_reduced_Data2S[-sid,]
dim(train)
dim(test)
head(test[1:12])
tail(test[1:12])
levels(test$Class.Label)
print(table(test$Class.Label))
print(prop.table(table(test$Class.Label)))

dim(train)
head(train[1:12])
tail(train[1:12])
levels(train$Class.Label)
print(table(train$Class.Label)) 
print(prop.table(table(train$Class.Label)))

###Balance the Dataset using SMOTE###
f_train <-filter(train, Class.Label == "alpha" | Class.Label == "PP")
levels(f_train$Class.Label) <- list(alpha="alpha", PP="PP")
f_train2 <-filter(train, Class.Label == "alpha" | Class.Label == "delta")
levels(f_train2$Class.Label) <- list(alpha="alpha", delta="delta")
f_train3 <-filter(train, Class.Label == "alpha" | Class.Label == "beta")
levels(f_train3$Class.Label) <- list(alpha="alpha", beta="beta")
dim(f_train)
dim(f_train2)
dim(f_train3)
head(f_train)
tail(f_train)
levels(f_train3$Class.Label)
as.factor(f_train3$Class.Label)
print(table(f_train2$Class.Label)) 
print(prop.table(table(f_train3$Class.Label)))

require(DMwR)
set.seed(1234)
###Synthetic Minority Over-Sampling Technique (SMOTE) Sampling
###This method is used to avoid overfitting when adding exact replicas of minority instances to the main dataset. For example, a subset of data from the minority class is taken. New synthetic similar instances are created and added to the original dataset. The count of each class records after applying                    
smote_train <- SMOTE(Class.Label ~., data=f_train,perc.over=800,perc.under=112.5)
print(table(smote_train$Class.Label))
print(prop.table(table(smote_train$Class.Label)))
dim(smote_train)
my_data1 = subset(smote_train, Class.Label == "PP")
my_data1$Class.Label


smote_train2 <- SMOTE(Class.Label ~., data=f_train2,perc.over=1595,perc.under=107)
print(table(smote_train2$Class.Label))
print(prop.table(table(smote_train2$Class.Label)))
dim(smote_train2)
my_data2 = subset(smote_train2, Class.Label == "delta")
my_data2$Class.Label

smote_train3 <- SMOTE(Class.Label ~., data=f_train3,perc.over=79,perc.under=228)
print(table(smote_train3$Class.Label)) 
print(prop.table(table(smote_train3$Class.Label)))
dim(smote_train3)
my_data3 = subset(smote_train3, Class.Label == "beta")
my_data3$Class.Label

my_data4 = subset(train, Class.Label == "alpha")
my_data4$Class.Label
BalanceTrain<-rbind(my_data4,my_data3,my_data2,my_data1)
dim(BalanceTrain)
levels(BalanceTrain$Class.Label)
BalanceTrain$Class.Label
print(table(BalanceTrain$Class.Label))
print(prop.table(table(BalanceTrain$Class.Label)))
write.csv(BalanceTrain, file = "BalanceTrain2.csv", row.names = FALSE)
BalanceTrain<-read.csv("BalanceTrain2.csv",header=T)



head(BalanceTrain)
tail(BalanceTrain)
head(test)
dim(BalanceTrain)
tail(test)
print(table(BalanceTrain$Class.Label))
print(prop.table(table(BalanceTrain$Class.Label)))

test.label<-test$Class.Label
###Building Models
require(randomForest)
require(ROCR)
model1 <- randomForest(Class.Label ~ ., data = BalanceTrain)
pred1 <- predict(model1, newdata = test)
cm<-table(pred1, test.label)   
(294 + 145 + 14+21) / nrow(test)   ###0.9978947
n = sum(cm) # number of instances
nc = nrow(cm) # number of classes
diag = diag(cm) # number of correctly classified instances per class 
rowsums = apply(cm, 1, sum) # number of instances per class
colsums = apply(cm, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1) 
"     precision    recall        f1
alpha 1.0000000 0.9966102 0.9983022
beta  1.0000000 1.0000000 1.0000000
delta 0.9333333 1.0000000 0.9655172
PP    1.0000000 1.0000000 1.0000000"

# Grid Search
seed<-1234
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(Class.Label ~ ., data=BalanceTrain, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
pred1.afterTune <- predict(rf_gridsearch, newdata = test)
table(pred1.afterTune, test.label)

confusionMatrix(data = pred1.afterTune, test.label)

(294 + 145 + 14+21) / nrow(test) ###0.9784615
table(test.label)
table(pred1.afterTune)
variable.imp<-varImp(rf_gridsearch)
variable.imp
tiff("Top10FeaturesImportanceDataset3",res=300,height=5,width=5,units="in")
plot(variable.imp, top = 10,ylab="Genes coding profile",main="Top 10 features")
dev.off()

require(dplyr)
require(dplyr)
train_new<-subset(T2D_reduced_Data2S, select=c("Class.Label","V5015","V1923","V4377" ,"V2673" ,"V4050","V5462","V1666" ,"V26473","V2503","V5317"))

#train_new <- train[,c("Class.Label","UBE2C.1","CCNB1.1","TOP2A ","HIST1H4C" ,"CDK1.1","RPL9.1","HMGB2.1" ,"HJURP","ZNF554","DNAJA1")]
train_new_<-train_new[,-1]
head(train_new)
write.csv(train_new, file = "GSE_data_10MostImpVarDataset3.csv", row.names = FALSE)
train_new<-read.csv("GSE_data_10MostImpVarDataset3.csv",header=T)
head(train_new)
dim(train_new)

###random forest on top10
train<-sample_frac(train_new, 0.70)
sid<-as.numeric(rownames(train)) # because rownames() returns character
test<-train_new[-sid,]
dim(train)
dim(test)
head(test[1:10])
tail(test[1:10])
levels(test$Class.Label)
print(table(test$Class.Label))
print(prop.table(table(test$Class.Label))) ###0.8513514


require(reshape)

train_new.m<-melt(train_new, id.var="Class.Label")

train_new.m
#****boxplot

ggplot(train_new.m, aes(x=variable,y=(value),fill=Class.Label))+ggtitle("Top 10 Variables' Class Distribution Dataset3")+theme_light()+geom_boxplot()

###*** Heatmap.2
require(gplots)
H <- train_new %>% group_by(Class.Label) %>% summarise_each(funs(mean))
H
H <- t(H)
H <- H[-1,]
H
H <- read.table("cellcycleheat3.txt")
H

heatmap.2(as.matrix(H),scale="row")
heatmap.2(as.matrix(H),scale="row",col=bluered)
heatmap.2(as.matrix(H),scale="row",col=bluered,Colv=T,tracecol=NA)###***

###***
pdf("outputHetmapDataset3.pdf",width=10,height=20)
heatmap.2(as.matrix(H),scale="row",col=bluered,Colv=FALSE,tracecol=NA)###***
dev.off()
test.label<-test$Class.Label
#### Model2: SVM
require("e1071")
model2 <- svm_model <- svm(Class.Label ~ ., data=BalanceTrain)
summary(model2)
pred2 <- predict(model2,newdata=test[,-96])
cm2<-table(pred2,test.label)
confusionMatrix(data = pred2, test.label)###Accuracy : 0.9298         


acc2<-(285 + 145 + 11+17) / nrow(test)
acc2   ###****0.9642105
print(table(test.label))
print(table(pred2))
n = sum(cm2) # number of instances
nc = nrow(cm2) # number of classes
diag = diag(cm2) # number of correctly classified instances per class 
rowsums = apply(cm2, 1, sum) # number of instances per class
colsums = apply(cm2, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"precision    recall        f1
alpha 0.9763314 0.9939759 0.9850746  #Balanced #0.9649123
beta  1.0000000 0.9042553 0.9497207
delta 0.7500000 1.0000000 0.8571429
PP    0.8421053 1.0000000 0.9142857"

" precision    recall        f1
alpha 0.9940828 0.9940828 0.9940828  #Imbalanced #acc=0.9298246
beta  1.0000000 0.8173077 0.8994709
delta 0.3333333 1.0000000 0.5000000
PP    0.4210526 1.0000000 0.5925926"
" precision    recall        f1
alpha 1.0000000 0.9941176 0.9970501 #Top10  ##Acc= 0.9929825
beta  1.0000000 0.9883721 0.9941520
delta 0.9166667 1.0000000 0.9565217
PP    0.9473684 1.0000000 0.9729730"

#### Model3:Naive Bayes:
model3 <- naiveBayes(Class.Label ~ ., data =BalanceTrain)
class(model3)
summary(model3)
print(model3)
pred3 <- predict(model3, newdata = test)
system.time(predict(model3,test))
cm3<-table(pred3,test.label)
confusionMatrix(data = pred3, test.label)
n = sum(cm3) # number of instances
nc = nrow(cm3) # number of classes
diag = diag(cm3) # number of correctly classified instances per class 
rowsums = apply(cm3, 1, sum) # number of instances per class
colsums = apply(cm3, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ####0.9452632
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"precision    recall        f1
alpha 0.9704142 0.9820359 0.9761905
beta  0.9882353 0.9767442 0.9824561
delta 0.7500000 0.8181818 0.7826087
PP    0.8421053 0.7619048 0.8000000"
print(table(test.label))
print(table(pred3))

control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"

nb_tune <- train(Class.Label~., data=BalanceTrain, method="nb",metric=metric, trControl=control)
nb_tune
pred3.afterTune <- predict(nb_tune, newdata = test)
system.time(predict(model3,test))
cm3<-table(pred3.afterTune,test.label)
n = sum(cm3) # number of instances
nc = nrow(cm3) # number of classes
diag = diag(cm3) # number of correctly classified instances per class 
rowsums = apply(cm3, 1, sum) # number of instances per class
colsums = apply(cm3, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ####***0.9852632
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"    precision    recall        f1
alpha 0.9829932 0.9965517 0.9897260
beta  0.9931034 1.0000000 0.9965398
delta 0.9333333 0.7368421 0.8235294
PP    1.0000000 0.9545455 0.9767442"
confusionMatrix(data = pred3.afterTune, test.label)

###_____________________________________________________________________________________________________________________________
#### Model4: K- Nearest Neighbors:
BalanceTrain[["Class.Lable"]] = factor(BalanceTrain[["Class.Label"]])
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(1234)
knn_fit <- train(Class.Lable ~., data = BalanceTrain, method = "knn",
                 trControl=trctrl,
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
knn_fit
plot(knn_fit)
head(test)
pred6_test <- predict(knn_fit, newdata = test)
pred6_test
system.time(predict(knn_fit,test))
cm4<-table(pred6_test,test.label)
n = sum(cm4) # number of instances
nc = nrow(cm4) # number of classes
diag = diag(cm4) # number of correctly classified instances per class 
rowsums = apply(cm4, 1, sum) # number of instances per class
colsums = apply(cm4, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ##########*******0.9810526
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"     precision    recall        f1
alpha 0.9704142 0.9820359 0.9761905
beta  0.9882353 0.9767442 0.9824561
delta 0.7500000 0.8181818 0.7826087
PP    0.8421053 0.7619048 0.8000000"
confusionMatrix(data = pred6_test, test.label)

####*****Decision Tree
control <- trainControl(method="repeatedcv", number=10, repeats=3,savePred =F)

set.seed(1234)
modelDT <- train(Class.Label~., data=BalanceTrain, method="rpart", trControl=control)
pred_dt <- predict(modelDT, test,type="raw")
pred_dt
system.time(predict(modelDT,test))
confusionMatrix(data = pred_dt, test.label)###Accuracy : 0.9473684               

cm5<-table(pred_dt,test.label)
n = sum(cm5) # number of instances
nc = nrow(cm5) # number of classes
diag = diag(cm5) # number of correctly classified instances per class 
rowsums = apply(cm5, 1, sum) # number of instances per class
colsums = apply(cm5, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ###0.9473684211
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)

# train the Model Averaged Neural Network
set.seed(1234)
modelANN <- train(Class.Label~., data=BalanceTrain, method="avNNet", trControl=control)
head(test)
dim(test)
pred_ANN <- predict(modelANN, test,type="raw")
pred_ANN
system.time(predict(modelANN,test))
confusionMatrix(data = pred_ANN, test.label)
cm6<-table(pred_ANN,test.label)
n = sum(cm6) # number of instances
nc = nrow(cm6) # number of classes
diag = diag(cm6) # number of correctly classified instances per class 
rowsums = apply(cm6, 1, sum) # number of instances per class
colsums = apply(cm6, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ###***0.9915789474
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"       precision    recall        f1
alpha 1.0000000 0.9941176 0.9970501
beta  1.0000000 1.0000000 1.0000000
delta 1.0000000 1.0000000 1.0000000
PP    0.9473684 1.0000000 0.9729730    "

"    precision    recall        f1
alpha 1.0000000 1.0000000 1.0000000   #Top10
beta  1.0000000 0.9883721 0.9941520
delta 0.9166667 1.0000000 0.9565217
PP    1.0000000 1.0000000 1.0000000"

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales) 
dat <- read.table("test_pre.txt",header=TRUE)
data.m <- melt(dat, id.vars=c('Class','Method'))
data.m           # sort columns


# melt the data frame for plotting

ggplot(data.m, aes(Method, value)) + geom_bar(aes(fill = Class), 
                                              width = 0.4, position = position_dodge(width=0.5), stat="identity") +  
  theme(legend.position="top", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank())+scale_fill_manual(values = c("#383838"	,"#666666","#949494","#B3B3B3"))+scale_y_continuous(labels=percent)+ggtitle("F1 Score in Percentage")+theme_light()     ####***Last Viz+xlab(label)





####RF caret
set.seed(1234)
modelRF <- train(Class.Label~., data=BalanceTrain, method="rf", trControl=control)

pred_RF <- predict(modelRF, test,type="raw")
system.time(predict(modelRF,test))
confusionMatrix(data = pred_RF, test.label)
cm7<-table(pred_RF,test.label)
n = sum(cm7) # number of instances
nc = nrow(cm7) # number of classes
diag = diag(cm7) # number of correctly classified instances per class 
rowsums = apply(cm7, 1, sum) # number of instances per class
colsums = apply(cm7, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ###0.9978947368
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"precision recall        f1
alpha 1.0000000   1.00 1.0000000 #0.9964912
beta  1.0000000   1.00 1.0000000
delta 0.9166667   1.00 0.9565217
PP    1.0000000   0.95 0.9743590"
##*** train the SVM model
set.seed(1234)
modelSVM <- train(Class.Label~., data=BalanceTrain, method="svmLinear", trControl=control)
pred_SVM <- predict(modelSVM, test,type="raw")
system.time(predict(modelSVM,test))
confusionMatrix(data = pred_SVM, test.label) ###Accuracy : 0.9768421
cm8<-table(pred_SVM,test.label)
n = sum(cm8) # number of instances
nc = nrow(cm8) # number of classes
diag = diag(cm8) # number of correctly classified instances per class 
rowsums = apply(cm8, 1, sum) # number of instances per class
colsums = apply(cm8, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"   precision       recall           f1
alpha 1.0000000000 0.9672131148 0.9833333333
beta  0.9929577465 0.9929577465 0.9929577465
delta 0.4615384615 1.0000000000 0.6315789474
PP    0.8800000000 1.0000000000 0.9361702128"
###** Train the NB model
set.seed(1234)
modelNB <- train(Class.Label~., data=BalanceTrain, method="nb", trControl=control)
pred_NB <- predict(modelNB, test,type="raw")
system.time(predict(modelNB,test))
confusionMatrix(data = pred_NB, test.label)
cm9<-table(pred_NB,test.label)
n = sum(cm9) # number of instances
nc = nrow(cm9) # number of classes
diag = diag(cm9) # number of correctly classified instances per class 
rowsums = apply(cm9, 1, sum) # number of instances per class
colsums = apply(cm9, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ###0.9642105263
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"         precision       recall           f1
alpha 0.9593220339 0.9964788732 0.9775474957
beta  0.9718309859 1.0000000000 0.9857142857
delta 0.9230769231 0.5217391304 0.6666666667
PP    1.0000000000 0.8333333333 0.9090909091"
###***Train the KNN model
set.seed(1234)
modelKNN <- train(Class.Label~., data=BalanceTrain, method="knn", trControl=control)
pred_KNN <- predict(modelKNN, test,type="raw")
system.time(predict(modelANN,test))
confusionMatrix(data = pred_KNN, test.label)
cm10<-table(pred_KNN,test.label)
n = sum(cm10) # number of instances
nc = nrow(cm10) # number of classes
diag = diag(cm10) # number of correctly classified instances per class 
rowsums = apply(cm10, 1, sum) # number of instances per class
colsums = apply(cm10, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy ###0.9789473684
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
KNN<-data.frame(precision, recall, f1)
"      precision       recall           f1
alpha 1.0000000000 0.9703947368 0.9849749583
beta  1.0000000000 1.0000000000 1.0000000000
delta 0.5384615385 0.8750000000 0.6666666667
PP    0.8400000000 1.0000000000 0.9130434783"
ibrary(ggplot2)
library(reshape2)
dat <- read.table("test_pre.txt",header=TRUE)
gg <- melt(dat,id=1:2)
ggplot(gg) +
  geom_bar(aes(x=Method, y=value, fill=Metric), stat="identity",
           position="dodge")+facet_wrap(~variable)+scale_fill_manual(values = c("#a1d99b", "#31a354"))



###scale_fill_manual(values = c("#006400","#a1d99b", "#31a354" ))
# collect resamples
results2 <- resamples(list(RF=modelRF, SVM=modelSVM, DT=modelDT,NB=modelNB,KNN=modelKNN,ANN=modelANN))
summary(results2)
# boxplots of results
