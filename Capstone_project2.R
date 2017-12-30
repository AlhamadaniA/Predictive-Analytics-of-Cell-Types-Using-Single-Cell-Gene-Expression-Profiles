devtools::install_github("hadley/readxl") # development version
library(readxl)
install.packages("rmarkdown")
require(rmarkdown)
require(ggpubr)
require(dplyr)
### Transposing and cleaning Non-diabetic dataset1
# Specify sheet with a number or name
NoNT2D_df<-read_excel("NoNT2D.xlsx", sheet = "Sheet3")
#read_excel("my-spreadsheet.xls", sheet = 2)
NoNT2D_df.T <- t(as.matrix(as.data.frame(NoNT2D_df)))
NoNT2D_df.T<-as.data.frame(NoNT2D_df.T)
head(NoNT2D_df.T[1:20], n=6)
tail(NoNT2D_df.T[1:20], n=6)
reponse<-read_excel("NoNT2D.xlsx", sheet = "Sheet6")
head(reponse)
write.csv(NoNT2D_df.T, file = "NoNT2D_df.T3.csv", row.names = FALSE)
NoNT2D_df.T3<-read.csv("NoNT2D_df.T3.csv",header=T)

NoNT2D_df.T3<-cbind(NoNT2D_df.T3,reponse)
head(NoNT2D_df.T3[39849:39852], n=6)
tail(NoNT2D_df.T3[39849:39852], n=6)

dim(NoNT2D_df.T3)
write.csv(NoNT2D_df.T3, file = "NoNT2D_df.T4.csv", row.names = FALSE)
NoNT2D_df.T3<-read.csv("NoNT2D_df.T4.csv",header=T)###*** Run this

###Step2:
### Transposing and cleaning diabetic dataset2
# Specify sheet with a number or name
T2D_df<-read_excel("T2D.xlsx", sheet = "Sheet1")
T2D_df.T <- t(as.matrix(as.data.frame(T2D_df)))
T2D_df.T<-as.data.frame(T2D_df.T)
head(T2D_df.T[1:20], n=6)
tail(T2D_df.T[1:20], n=6)
reponse2<-read_excel("T2D.xlsx", sheet = "Sheet2")
head(reponse2)
write.csv(T2D_df.T, file = "T2D_df.T1.csv", row.names = FALSE)
T2D_df.T1<-read.csv("T2D_df.T1.csv",header=T)

T2D_df.T1<-cbind(T2D_df.T1,reponse2)
dim(T2D_df.T1)
T2D_df.T1[,39854]<-NULL
head(T2D_df.T1[39849:39853], n=6)
tail(T2D_df.T1[39849:39853], n=6)

write.csv(T2D_df.T1, file = "T2D_df.T3.csv", row.names = FALSE)
T2D_df.T1<-read.csv("T2D_df.T3.csv",header=T)###*** Run This
###===========================================================================================
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
str(NoNT2D_df.T3)
dim(NoNT2D_df.T3)
head(NoNT2D_df.T3$Cell.Type)
names(NoNT2D_df.T3)
class(NoNT2D_df.T3)
NoNT2D_df.T3$Cell.Type
is.factor(NoNT2D_df.T3$Cell.Type)
###Convert the response from charachter to factor by index number
NoNT2D_df.T3[,39852] <- as.factor(NoNT2D_df.T3[,39852])
###Verify and check levels or factor types
is.factor(NoNT2D_df.T3[,39852])
levels(NoNT2D_df.T3$Cell.Type)
target<-NoNT2D_df.T3$Cell.Type
###Further explore the data
glimpse(NoNT2D_df.T3[1:20])
###***Check the balance of the class***
print(table(NoNT2D_df.T3$Cell.Type))
print(prop.table(table(NoNT2D_df.T3$Cell.Type))) ###There is class imbalance***
### Filter Data: Remove genes with low expression
NoNT2D_df.T3$Cell.Type<-NULL

myfunction <- function(x){
  
  {    P<-mean(x)
  }
  return(P)
}
alist=list() 
blist=list() 
for(i in 1:length(NoNT2D_df.T3)){
P<-myfunction(NoNT2D_df.T3[,i])

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

NoNT2D_df.T3Subset<- subset(NoNT2D_df.T3, select = desiredVector)
dim(NoNT2D_df.T3Subset)   ####24744
head(NoNT2D_df.T3Subset[1:12])
tail(NoNT2D_df.T3Subset[1:12])
head(NoNT2D_df.T3Subset[1:4])

dfNew<-cbind(target,NoNT2D_df.T3Subset)
dim(dfNew)
head(dfNew[1:12])
tail(dfNew[1:12])
levels(dfNew$target)
print(table(dfNew$target))
write.csv(dfNew, file = "dfNew.csv", row.names = FALSE)
dfNew<-read.csv("dfNew.csv",header=T)
dim(dfNew)
#dfNew.scaled <- as.data.frame(scale(dfNew[,-1]))
#str(dfNew.scaled)
#tail(dfNew.scaled[1:12])
sink("AnovaTest3.doc")
for(j in 2:ncol(dfNew)){
  #### create side by side boxplot of the data:
  boxplot(dfNew[,j]~dfNew$target)
  ###Create an ANOVA table:
  
  columns <- names(dfNew[j])
  anovaresult<- anova(aov(dfNew[,j]~target,data=dfNew)) 
  ###Conduct a Tukey's multiple comparison precedure:
  posthocresult <-TukeyHSD(aov(dfNew[,j]~target,data=dfNew))
  if(anovaresult$Pr[1] < 0.05){ 
    print(anovaresult) 
    flist = list(flist, j)}
  print(columns)
  #print(posthocresult)
}    
sink()

length(unlist(flist))
desiredVector2<-unlist(flist, use.names=FALSE)
NoNT2D_df.T3Subset2<- subset(dfNew, select = desiredVector2)
NoNT2D_Anova<-cbind(target,NoNT2D_df.T3Subset2)
dim(NoNT2D_Anova)
head(NoNT2D_Anova[1:12])
tail(NoNT2D_Anova[1:12])
levels(NoNT2D_Anova$target)
print(table(NoNT2D_Anova$target))
###Removing replicated variables:
set.seed(1234)
dim(NoNT2D_Anova)
df2 = cor(NoNT2D_Anova[,-1])
hc = findCorrelation(df2, cutoff=0.99) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = NoNT2D_Anova[,-1][,-c(hc)]
dim(reduced_Data)
head(reduced_Data[,1:12])
colnames(reduced_Data)
NoNT2D_reduced_Data<-cbind(target,reduced_Data)
dim(NoNT2D_reduced_Data) ###***6521
head(NoNT2D_reduced_Data[,1:12])
write.csv(NoNT2D_reduced_Data, file = "NoNT2D_reduced_Data.csv", row.names = FALSE)
NoNT2D_reduced_Data<-read.csv("NoNT2D_reduced_Data.csv",header=T)


#Lasso Regression Step
set.seed(1234)
cv.glmnet.fit <- cv.glmnet(as.matrix(NoNT2D_reduced_Data[,-1]),NoNT2D_reduced_Data[,1],alpha = 0.5, type.multinomial ="grouped",type.measure = "class",nlambda =200, family ="multinomial" )
plot(cv.glmnet.fit)
tmp<-coef(cv.glmnet.fit,s="lambda.min")
T<-as.matrix(tmp[[1]])
coefs<-data.frame(T[(abs(T[,1]) >0),])
features<-row.names(coefs)
features
idx<-match(features,colnames(NoNT2D_reduced_Data),nomatch = 0)
NoNT2D_reduced_Data2<-NoNT2D_reduced_Data[,idx]
NoNT2D_reduced_Data2$target<-NoNT2D_reduced_Data$target
head(NoNT2D_reduced_Data2[125:135])
tail(NoNT2D_reduced_Data2[125:135])
dim(NoNT2D_reduced_Data2)
colnames(NoNT2D_reduced_Data2)[135] <- "Class.Label"
dput(names(NoNT2D_reduced_Data2))
colnames(NoNT2D_reduced_Data2)<-c("V5775", "V6127", "V92", "V152", "V170", "V326", "V386", 
                                  "V437", "V440", "V457", "V523", "V861", "V884", "V914", 
                                  "V942", "V980", "V1163", "V1666", "V1921", "V1923", 
                                  "V2001", "V2118", "V2167", "V2252", "V2297", "V2475", 
                                  "V2503", "V2514", "V2673", "V2685", "V2740", "V3026", 
                                  "V3057", "V3117", "V3260", "V3269", "V3274", "V3282", 
                                  "V3304", "V3415", "V3545", "V3712", "V3714", "V3784", 
                                  "V3826", "V4015", "V4050", "V4070", "V4207", "V4377", 
                                  "V4415", "V4416", "V4494", "V4548", "V4553", "V4580", 
                                  "V4597", "V4824", "V4960", "V5015", "V5254", "V5292", 
                                  "V5316", "V5317", "V5340", "V5448", "V5462", "V5519", 
                                  "V5804", "V6344", "V6525", "V6758", "V6958", "V7068", 
                                  "V7176", "V7363", "V7432", "V7435", "V7447", "V7496", 
                                  "V8205", "V8223", "V8282", "V8747", "V9100", "V9241", 
                                  "V9418", "V9513", "V9530", "V9947", "V10489", "V11259", 
                                  "V11516", "V11541", "V11922", "V12156", "V12163", "V12180", 
                                  "V12264", "V12636", "V12658", "V12808", "V12904", "V12916", 
                                  "V13063", "V13187", "V13252", "V13656", "V13945", "V14315", 
                                  "V14555", "V14806", "V14941", "V15145", "V15856", "V17262", "V17828", 
                                  "V17914", "V17954", "V18060", "V18200", "V19164", "V21358", "V21602", 
                                  "V21646", "V21852", "V22129", "V22869", "V26473", "V26912", "V26951", 
                                  "V27007", "V27723", "V28066", "Class.Label")
write.csv(NoNT2D_reduced_Data2, file = "NoNT2D_reduced_Data2.csv", row.names = FALSE)
NoNT2D_reduced_Data2<-read.csv("NoNT2D_reduced_Data2.csv",header=T)
dim(NoNT2D_reduced_Data2)
ad<-NoNT2D_reduced_Data2 %>% count(Label = factor(Class.Label), Class.Label = factor(Class.Label)) %>% 
  ungroup() %>%
  mutate(percentage = round(prop.table(n) * 100),digits=2) %>% 
  ggplot(aes(x = Class.Label, y = percentage, fill = Class.Label)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  ylim(0, 100)+
  geom_text(aes(y = percentage + .5,    # nudge above top of bar
                label = paste0(percentage, '%'),size=14),    # prettify
            position = position_dodge(width = .9), 
            size = 5)+theme_light() +ggtitle("Imbalanced Class 'Endocrine Pancreas: Islets of Langerhans' Dataset2")
ad+theme(axis.text=element_text(size=12),
           axis.title=element_text(size=14,face="bold"))
###*** Training and Testing***###
dim(NoNT2D_reduced_Data2)


#colnames(NoNT2D_reduced_Data2)[135] <- "Class.Label"

###*** Training and Testing***###
dim(NoNT2D_reduced_Data2)
NoNT2D_reduced_Data2S<-scale(NoNT2D_reduced_Data2[,-135])
NoNT2D_reduced_Data2S<-as.data.frame(NoNT2D_reduced_Data2S)
NoNT2D_reduced_Data2S<-cbind(NoNT2D_reduced_Data2S,NoNT2D_reduced_Data2$Class.Label)
colnames(NoNT2D_reduced_Data2S)[135] <- "Class.Label"
colnames(NoNT2D_reduced_Data2S)
class(NoNT2D_reduced_Data2S)

train<-sample_frac(NoNT2D_reduced_Data2S, 0.70)
sid<-as.numeric(rownames(train)) # because rownames() returns character
test<-NoNT2D_reduced_Data2S[-sid,]
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

install.packages("DMwR")
require(DMwR)
set.seed(1234)
###Synthetic Minority Over-Sampling Technique (SMOTE) Sampling
###This method is used to avoid overfitting when adding exact replicas of minority instances to the main dataset. For example, a subset of data from the minority class is taken. New synthetic similar instances are created and added to the original dataset. The count of each class records after applying                    
smote_train <- SMOTE(Class.Label ~., data=f_train,perc.over=1099,perc.under=110)
print(table(smote_train$Class.Label))
print(prop.table(table(smote_train$Class.Label)))
dim(smote_train)
my_data1 = subset(smote_train, Class.Label == "PP")
my_data1$Class.Label


smote_train2 <- SMOTE(Class.Label ~., data=f_train2,perc.over=1300,perc.under=108)
print(table(smote_train2$Class.Label))
print(prop.table(table(smote_train2$Class.Label)))
dim(smote_train2)
my_data2 = subset(smote_train2, Class.Label == "delta")
my_data2$Class.Label

smote_train3 <- SMOTE(Class.Label ~., data=f_train3,perc.over=90,perc.under=213)
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
write.csv(BalanceTrain, file = "BalanceTrain.csv", row.names = FALSE)
BalanceTrain<-read.csv("BalanceTrain.csv",header=T)


head(BalanceTrain)
tail(BalanceTrain)
head(test)
tail(test)
print(table(test$Class.Label))
print(prop.table(table(test$Class.Label)))

test.label<-test$Class.Label
test$Class.Label<-NULL
###Building Models
require(randomForest)
require(ROCR)
model1 <- randomForest(Class.Label ~ ., data = BalanceTrain)
pred1 <- predict(model1, newdata = test)
cm<-table(pred1, test.label)   
(175 + 110 + 17+20) / nrow(test)   ###0.9938462
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


# Grid Search
seed<-1234
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(Class.Label ~ ., data=train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
pred1.afterTune <- predict(rf_gridsearch, newdata = test)
table(pred1.afterTune, test.label)

confusionMatrix(data = pred1.afterTune, test.label)

(188 + 111 + 12+12) / nrow(test) ###0.9938462
table(test.label)
cm<-table(pred1.afterTune,test.label)
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


variable.imp<-varImp(rf_gridsearch)
variable.imp
tiff("Top10FeaturesImportanceDataset2",res=300,height=5,width=5,units="in")
plot(variable.imp, top = 10,ylab="Genes coding profile",main="Top 10 features")
dev.off()

train_new<-subset(NoNT2D_reduced_Data2S, select=c("Class.Label","V1923","V5462","V2673" ,"V5317" ,"V2503","V2252","V1921" ,"V26473","V1666","V3784"))

write.csv(train_new, file = "GSE_data_10MostImpVarDataset2.csv", row.names = FALSE)
train_new<-read.csv("GSE_data_10MostImpVarDataset2.csv",header=T)
#subset(train, select=c("Class.Label","UBE2C.1","CCNB1.1","TOP2A ","HIST1H4C" ,"CDK1.1","RPL9.1","HMGB2.1" ,"HJURP","ZNF554","DNAJA1"))
#train_new <- train[,c("Class.Label","UBE2C.1","CCNB1.1","TOP2A ","HIST1H4C" ,"CDK1.1","RPL9.1","HMGB2.1" ,"HJURP","ZNF554","DNAJA1")]
train_new_<-train_new[,-1]
head(train_new)
write.csv(train_new, file = "GSE_data_10MostImpVarDataset2.csv", row.names = FALSE)
train_new<-read.csv("GSE_data_10MostImpVarDataset2.csv",header=T)
head(train_new)
dim(train_new)
require(reshape)

train_new.m<-melt(train_new, id.var="Class.Label")

train_new.m
#****boxplot

ggplot(train_new.m, aes(x=variable,y=(value),fill=Class.Label))+ggtitle("Top 10 Variables' Class Distribution Dataset2")+theme_light()+geom_boxplot()

###*** Heatmap.2
require(gplots)
H <- train_new %>% group_by(Class.Label) %>% summarise_each(funs(mean))
H
H <- t(H)
H <- H[-1,]
H
H <- read.table("cellcycleheat2.txt")
H

heatmap.2(as.matrix(H),scale="row")
heatmap.2(as.matrix(H),scale="row",col=bluered,Colv=T,tracecol=NA)

###***
pdf("outputHetmapDataset2.pdf",width=10,height=20)
heatmap.2(as.matrix(H),scale="row",col=bluered,Colv=FALSE,tracecol=NA)###***
dev.off()

#### Model2: SVM
require("e1071")
model2 <- svm_model <- svm(Class.Label ~ ., data=train_new)
summary(model2)
pred2 <- predict(model2,newdata=test[,-135])
cm2<-table(pred2,test.label)
confusionMatrix(data = pred2, test.label)

acc2<-(183 + 105 + 4+11) / nrow(test)
acc2
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
"Balanced precision    recall        f1
alpha    0.9824561 1.0000000 0.9911504
beta     1.0000000 0.9285714 0.9629630
delta    0.8888889 1.0000000 0.9411765
PP       0.7142857 1.0000000 0.8333333"

"Imbalanced precision recall        f1
alpha       0.9736842 1.0000 0.9866667
beta        0.9692308 1.0000 0.9843750
delta       0.5555556 1.0000 0.7142857
PP          1.0000000 0.4375 0.6086957"

"Top10 precision    recall        f1
alpha 1.0000000 0.9827586 0.9913043
beta  0.9846154 1.0000000 0.9922481
delta 0.7777778 1.0000000 0.8750000
PP    0.8571429 0.7500000 0.8000000"

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
accuracy ###0.9641026
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
" precision    recall        f1       ###balanced
alpha 0.9561404 0.9909091 0.9732143
beta  0.9846154 0.9846154 0.9846154
delta 1.0000000 1.0000000 1.0000000
PP    0.8571429 0.5454545 0.6666667"
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
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
confusionMatrix(data = pred3.afterTune, test.label)###0.9795

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
accuracy ##########*******0.9846154
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
confusionMatrix(data = pred6_test, test.label)
"  precision    recall        f1
alpha 0.9897959 0.9965753 0.9931741
beta  1.0000000 0.9731544 0.9863946
delta 0.7333333 0.7333333 0.7333333
PP    0.9047619 1.0000000 0.9500000"
####*****Decision Tree
control <- trainControl(method="repeatedcv", number=10, repeats=3,savePred =F)

set.seed(1234)
modelDT <- train(Class.Label~., data=BalanceTrain, method="rpart", trControl=control)
pred_dt <- predict(modelDT, test,type="raw")
pred_dt
system.time(predict(modelDT,test))
confusionMatrix(data = pred_dt, test.label)###   Accuracy : 0.9476923  
cm5<-table(pred_dt,test.label)
n = sum(cm5) # number of instances
nc = nrow(cm5) # number of classes
diag = diag(cm5) # number of correctly classified instances per class 
rowsums = apply(cm5, 1, sum) # number of instances per class
colsums = apply(cm5, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes
accuracy = sum(diag) / n 
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
" precision recall   f1
alpha         1 1.0000 1.00
beta          1 1.0000 1.00
delta         1 0.5625 0.72
PP            0    NaN  NaN"
# train the Model Averaged Neural Network
set.seed(1234)
modelANN <- train(Class.Label~., data=BalanceTrain, method="avNNet", trControl=control)
head(test)
dim(test)
test<-test[,-135]
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
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"    precision recall f1
alpha         1      1  1
beta          1      1  1
delta         1      1  1
PP            1      1  1"
###measurement plot
dat <- read.table("test_pre.txt",header=TRUE)
head()

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




# plot everything
ggplot(data.m, aes(Names, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity")
gg <- melt(dat,id=1:2)
ggplot(gg) +
  geom_bar(aes(x=Method, y=F1Score, fill=variable), stat="identity",
           position="dodge")+facet_wrap(~variable)+scale_fill_manual(values = c("#a1d99b", "#31a354"))




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
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
"      precision    recall        f1
alpha 1.0000000 0.9913043 0.9956332
beta  1.0000000 1.0000000 1.0000000
delta 1.0000000 1.0000000 1.0000000
PP    0.8571429 1.0000000 0.9230769"
##*** train the SVM model
set.seed(1234)
modelSVM <- train(Class.Label~., data=BalanceTrain, method="svmLinear", trControl=control)
pred_SVM <- predict(modelSVM, test,type="raw")
system.time(predict(modelSVM,test))
confusionMatrix(data = pred_SVM, test.label)
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
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)

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
accuracy 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1)
# collect resamples
results2 <- resamples(list(RF=modelRF, SVM=modelSVM, DT=modelDT,NB=modelNB,KNN=modelKNN,ANN=modelANN))
summary(results2)
# boxplots of results
bwplot(results2,Metric="roc")

