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
set.seed(1234)

GSE_data<-read.csv("GSE64016_H1andFUCCI_normalized_EC.csv",header=T)
str(GSE_data)
dim(GSE_data)
names(GSE_data)
class(GSE_data)
levels(GSE_data$Class.Label)
head(GSE_data[1:20], n=4)
tail(GSE_data[1:15], n=4)
glimpse(GSE_data[1:20])
print(table(GSE_data$Class.Label))
print(prop.table(table(GSE_data$Class.Label)))
hist(GSE_data$SGMS2)
GSE_data$Sample.Name <-NULL
target<-GSE_data$Class.Label
#GSE_data1<-GSE_data[-1]
View(GSE_data)
#dim(GSE_data1)
#GSE_data$Class.Label<-NULL
str(GSE_data)
# One Way Anova (Completely Randomized Design)
#fit_aov <- aov(Class.Label ~., data=GSE_data)
#summary(fit_aov)

####***Remove genes with low expressions
myfunction <- function(x,y){
  
  {    P<-(x*100)/y
}
  return(P)
}
alist=list() 
blist=list() 
clist=list() 
dlist=list() 
elist=list() 
flist=list()
glist=list()
for(i in 2:length(GSE_data)){
P_G1 <-myfunction(nrow(subset(GSE_data[1:92,],GSE_data[1:92,][,i]!="0")),91)
P_S<-myfunction(nrow(subset(GSE_data[93:172,], GSE_data[93:172,][,i]!="0")),80)
P_G2<-myfunction(nrow(subset(GSE_data[173:248,], GSE_data[173:248,][,i]!="0")),76)
P_0<-nrow(subset(GSE_data[1:248,], GSE_data[1:248,][,i]!=0))
if (P_G1>50) {
  print(paste("The column is", i, "keep G1> 50%" ))
  alist = list(alist, i) # combining them in a list
  
  #result <- do.call("cbind", extracted_cols) # trying to cbind the third colum
  
}else if(P_S>50)
  {
  print(paste("The column is", i, "keep S > 50%" ))
  blist = list(blist, i) # combining them in a list
  
  #result2<- do.call("cbind", GSE_data[i])
  }else if(P_G2>50)
{
  print(paste("The column is", i, "keep G2 > 50%" ))
  #result3<- do.call("cbind", GSE_data[i])
    clist = list(clist, i) # combining them in a list
    
}

 else if(P_0==0)
 {
   print(paste("The column is", i, "all 0's skip and omit" ))
   dlist = list(dlist, i) # combining them in a list
   
 }else 
   print(paste("The column is", i, "skip and omit" ))
elist = list(elist, i) # combining them in a list

}
desiredList<-do.call(c, list(alist, blist,clist))
undesiredList<-do.call(c, list(dlist, elist))

desiredVector<-unlist(desiredList, use.names=FALSE)
undesiredVector<-unlist(undesiredList, use.names=FALSE)

head(desiredVector)
tail(desiredVector)
head(undesiredVector)
tail(undesiredVector)
length(desiredVector)
length(undesiredVector)

GSE_dataSubset<- subset(GSE_data, select = desiredVector)
dim(GSE_dataSubset)
head(GSE_dataSubset[1:12])
tail(GSE_dataSubset[1:12])
head(GSE_data[1:4])

df22<-cbind(GSE_data[1],GSE_dataSubset)
dim(df22)
head(df22[1:12])
tail(df22[1:12])
levels(df22$Class.Label)
print(table(df22$Class.Label))
write.csv(df22, file = "GSE_data22.csv", row.names = FALSE)
df22<-read.csv("GSE_data22.csv",header=T)
dim(df22)
# One Way Anova (Completely Randomized Design)
sink("AnovaTest.doc")
for(j in 2:ncol(df22)){
#### create side by side boxplot of the data:
boxplot(df22[,j]~df22$Class.Label)
###Create an ANOVA table:

columns <- names(df22[j])
anovaresult<- anova(aov(df22[,j]~Class.Label,data=df22)) 
###Conduct a Tukey's multiple comparison precedure:
posthocresult <-TukeyHSD(aov(df22[,j]~Class.Label,data=df22))
if(anovaresult$Pr[1] < 0.05){ 
print(anovaresult) 
flist = list(flist, j)}
print(columns)
#print(posthocresult)
}    
sink()

length(unlist(flist))
desiredVector2<-unlist(flist, use.names=FALSE)
GSE_dataSubset2<- subset(df22, select = desiredVector2)
df222<-cbind(df22[1],GSE_dataSubset2)
dim(df222)
head(df222[1:12])
tail(df222[1:12])
levels(df222$Class.Label)
print(table(df222$Class.Label))
write.csv(df222, file = "GSE_data222.csv", row.names = FALSE)
df222<-read.csv("GSE_data222.csv",header=T)
###Removing replicated variables:
set.seed(1234)
dim(df222)
target<-df222[,1]
df2 = cor(df222[,-1])
cor(df222$EIF3J,df222$EIF3J.1)
hc = findCorrelation(df2, cutoff=0.99) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = df222[,-1][,-c(hc)]
dim(reduced_Data)
head(reduced_Data[,1:12])
colnames(reduced_Data)
reduced_Data<-cbind(target,reduced_Data)
colnames(reduced_Data)[1] <- "Class.Label"
write.csv(reduced_Data, file = "GSE_reduced_Anova.csv", row.names = FALSE)
reduced_Data<-read.csv("GSE_reduced_Anova.csv",header=T)
dim(reduced_Data)
dput(names(reduced_Data))
colnames(train_new)<-c("Class.Label", "UBE2C", "CCNB1", "TOP2A", "HIST1H4C", "CDK1", 
                       "RPL9", "HMGB2", "HJURP", "ZNF554", "DNAJA1")
#lasso
set.seed(1234)
x<-as.matrix(reduced_Data[,-1])
y<-reduced_Data[,1]
cv.glmnet.fit <- cv.glmnet(x,y,alpha = 0.5, type.multinomial ="grouped",type.measure = "class",nlambda =200, family ="multinomial" )
plot(cv.glmnet.fit)
tmp<-coef(cv.glmnet.fit,s="lambda.min")
T<-as.matrix(tmp[[1]])
coefs<-data.frame(T[(abs(T[,1]) >0),])
features<-row.names(coefs)
features
idx<-match(features,colnames(df222),nomatch = 0)
df222_GSE<-df222[,idx]
df222_GSE$Class.Label<-df222$Class.Label
dim(df222_GSE)
head(df222_GSE[247:279])
write.csv(df222_GSE, file = "reducedData_logistic.csv", row.names = FALSE)
df222_GSE<-read.csv("reducedData_logistic.csv",header=T)
dput(names(df222_GSE))

opt_lambda <- cv.glmnet.fit$lambda.min
opt_lambda
fit <- cv.glmnet.fit$glmnet.fit
summary(fit)
y_predicted <- predict(fit, s = opt_lambda, newx = as.matrix(reduced_Data[,-1]))


###Removing replicated variables:
set.seed(1234)
dim(df222_GSE)

barggplot_class<-ggplot(df222_GSE, aes(x = factor(Class.Label),fill=Class.Label)) + geom_bar(stat = "count")
print(barggplot_class + ggtitle("Balanced Class Label Distribution of Dataset1"))

df222_GSE %>% count(Label = factor(Class.Label), Class.Label = factor(Class.Label)) %>% 
  ungroup() %>%
  mutate(percentage = round(prop.table(n) * 100),digits=2) %>% 
  ggplot(aes(x = Class.Label, y = percentage, fill = Class.Label)) +
  geom_bar(stat = 'identity', position = 'dodge')  +scale_fill_brewer(palette="Dark2")+
  ylim(0, 100)+
  geom_text(aes(y = percentage + .5,    # nudge above top of bar
                label = paste0(percentage, '%')),    # prettify
            position = position_dodge(width = .9), 
            size = 3)+ggtitle("Balanced Class 'Cell Cycle Types' ")

aad<-df222_GSE %>% count(Label = factor(Class.Label), Class.Label = factor(Class.Label)) %>% 
  ungroup() %>%
  mutate(percentage = round(prop.table(n) * 100),digits=2) %>% 
  ggplot(aes(x = Class.Label, y = percentage, fill = Class.Label)) + 
  geom_bar(stat = 'identity', position = 'dodge') +scale_fill_brewer(palette="Dark2")+
  ylim(0, 100)+
  geom_text(aes(y = percentage + .5,    # nudge above top of bar
                label = paste0(percentage, '%'),size=14),    # prettify
            position = position_dodge(width = .9), 
            size = 5)+theme_light() +ggtitle("Balanced Class Embryonic Cell Types Data")
aad+theme(axis.text=element_text(size=12),
           axis.title=element_text(size=14,face="bold"))

library(ggplot2)
library(scales)
ggplot(mtcars, aes(x = factor(hp))) +  
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels = percent, limits = c(0,1))

###removing the intercept column), store the independent variable as y, and create a vector of lambda values.
"x <- model.matrix(Class.Label~., df222)[,-1]
x<-model.matrix(Class.Label~.-1,data=df222) 
y <- unclass(df222$Class.Label) %>% as.numeric 
y <- df222$Class.Label

lambda <- 10^seq(10, -2, length = 100)
head(x[1:12])
levels(y)
# fit model
fit.lasso=glmnet(x,y,alpha=1)
summary(fit.lasso)
plot(fit.lasso,xvar="lambda",label=TRUE)
coef(fit.lasso, s = 0.3) # s is the value of lambda
crossval <-  cv.glmnet(x = x, y = y)
plot(crossval)
penalty <- crossval$lambda.min #optimal lambda
penalty #minimal shrinkage
fit1 <-glmnet(x = x, y = y, alpha = 1, lambda = penalty ) #estimate the model with that
coef12<-coef(fit1)



plot(fit.lasso,xvar="dev",label=TRUE)
cv.lasso=cv.glmnet(x,y)
plot(cv.lasso)
coef(cv.lasso)
fit.lasso=glmnet(x,y,alpha=1)
names(fit.lasso)
fit.lasso
fit.lasso$df
#for(i in 2:length(GSE_data)){
output<-data.frame(matrix(nrow=4928, ncol=4))
names(output)=c("Estimate", " Std. Error", " z value", " Pr(>|z|)")
for (i in 2:length(df222)){
  resultLasso<-(glm(df222$Class.Label ~ as.matrix(df222[,i]) , family=binomial("logit"),maxit = 100))
  row.names(output)<-row.names(4928)

  output[i,1]<-coef(summary(resultLasso))[2,1]
  output[i,2]<-coef(summary(resultLasso))[2,2]
  output[i,3]<-coef(summary(resultLasso))[2,3]
  output[i,4]<-coef(summary(resultLasso))[2,4] 
  if (output[i,4]<0.05){
    glist=list(glist,i)
  }
}
length(unlist(glist))
desiredVector3<-unlist(glist, use.names=FALSE)
GSE_dataSubset3<- subset(df222, select = desiredVector3)
df2223<-cbind(df222[1],GSE_dataSubset3)
dim(df2223)
head(df2223[1:12])
tail(df2223[1:12])
levels(df2223$Class.Label)
print(table(df2223$Class.Label))
write.csv(df2223, file = "GSE_data2223.csv", row.names = FALSE)
df2223<-read.csv("GSE_data2223.csv",header=T)
getwd()
df2223$Class.Label<-as.numeric(df2223$Class.Label)
df2223$Class.Label
df2223$Class.Label<-as.factor(df2223$Class.Label)
"
###**********Vizualization **************************************************************************************************************
options(scipen=999)  # turn-off scientific notation like 1e+48
theme_set(theme_bw())  # pre-set the bw theme.
gg <- ggplot(df222_GSE, aes(x= MCM3, y=MCM6)) + 
  geom_point(aes(col=Class.Label, size=MCM6)) + 
  geom_smooth(method="loess", se=F) + 
 # xlim(c(0, 0.1)) + 
 # ylim(c(0, 500000)) + 
  labs(subtitle="MCM3 Vs MCM6", 
       y="Class.Label", 
       x="MCM3", 
       title="Scatterplot", 
       caption = "Source: GSE_data2223")

plot(gg)


Class.Label_select <- df222_GSE[df222_GSE$Class.Label %in% c("G1", "S", "G2"), ]

# Scatterplot
theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(Class.Label_select, aes(MCM3, MCM6)) + 
  labs(subtitle="Gene Classification: MCM3 vs MCM6",
       title="Bubble chart")

g + geom_jitter(aes(col=Class.Label, size=MCM3)) + 
  geom_smooth(aes(col=Class.Label), method="lm", se=F)


###****************End of Vizualization **************************************************************************
###6. Training and Testing
summary(df222_GSE)
set.seed(1234)
scaled_df222_GSE<-scale(df222_GSE[,-279])
scaled_GSE<-as.data.frame(scaled_df222_GSE)
scaled_GSE<-cbind(scaled_GSE,df222_GSE$Class.Label)
colnames(scaled_GSE)[279] <- "Class.Label"
colnames(scaled_GSE)
class(scaled_GSE)
dim(scaled_GSE)
write.csv(scaled_GSE, file = "scaled_GSE.csv", row.names = FALSE)
scaled_GSE<-read.csv("scaled_GSE.csv",header=T)

library(dplyr)
train<-sample_frac(scaled_GSE, 0.70)
sid<-as.numeric(rownames(train)) # because rownames() returns character
test<-scaled_GSE[-sid,]
dim(train)
dim(test)
head(test[1:12])
tail(test[1:12])
levels(test$Class.Label)
print(table(test$Class.Label))
dim(train)
head(train[1:12])
tail(train[1:12])
levels(train$Class.Label)
print(table(train$Class.Label))
###******************************************************************************
###Building the models
require(randomForest)
test.label<-test[,279]
length(test.label)
#test<-test[,-80]
dim(test)
train.label<-train[,279]
length(train.label)

model1 <- randomForest(Class.Label ~ ., data = train)
pred1 <- predict(model1, newdata = test)
table(pred1, test.label)
####test the accuracy 
(26 + 21 + 18) / nrow(test)###*** (26 + 21 + 18) / nrow(test) [1] 0.8783784***
print(table(test.label))
print(table(pred1))

###Extract feature (X) and target (y) columns
x<-train[,-279]
y<-train[,279]  
dim(x)
length(y)

### Create trainControl object(with 10 fold Cross Validation): myControl
# Create model with default paramters
library(caret)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
seed <- 1234
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class.Label~., data= test, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)   "Accuracy   Kappa    
                    0.8314153  0.7414755"

# Grid Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(Class.Label ~ ., data=train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
pred1.afterTune <- predict(rf_gridsearch, newdata = test)
table(pred1.afterTune, test.label)

confusionMatrix(data = pred1.afterTune, test.label)

(27 + 21 + 20) / nrow(test) ###0.9189189
table(test.label)
table(pred1.afterTune)
variable.imp<-varImp(rf_gridsearch)
variable.imp
tiff("Top10FeaturesImportanceScaled",res=300,height=5,width=5,units="in")

plot(variable.imp, top = 10,ylab="Genes coding profiles",main="Top 10 features")
dev.off()

str(scaled_GSE)
require(dplyr)
train_new<-subset(scaled_GSE, select=c("Class.Label","UBE2C","HIST1H4C","CDK1" ,"SRSF7" ,"CCNB1" ,"TOP2A","HJURP","BUB3" ,"PMAIP1","TLK1"))
train<-sample_frac(train_new, 0.70)
sid<-as.numeric(rownames(train)) # because rownames() returns character
test<-train_new[-sid,]
dim(train)
dim(test)
head(test[1:10])
tail(test[1:10])
levels(test$Class.Label)
print(table(test$Class.Label))
dim(train)
head(train[1:10])
tail(train[1:10])
levels(train$Class.Label)
print(table(train$Class.Label))
###Building the models
require(randomForest)
test.label<-test[,279]
length(test.label)
#test<-test[,-80]
dim(test)
train.label<-train[,279]
length(train.label)

model1 <- randomForest(Class.Label ~ ., data = train)
pred1 <- predict(model1, newdata = test)
table(pred1, test.label)
####test the accuracy 
(22 + 18 + 23) / nrow(test)###*** (26 + 21 + 18) / nrow(test) [1] 0.8783784***
print(table(test.label))
print(table(pred1))



#subset(train, select=c("Class.Label","UBE2C.1","CCNB1.1","TOP2A ","HIST1H4C" ,"CDK1.1","RPL9.1","HMGB2.1" ,"HJURP","ZNF554","DNAJA1"))
#train_new <- train[,c("Class.Label","UBE2C.1","CCNB1.1","TOP2A ","HIST1H4C" ,"CDK1.1","RPL9.1","HMGB2.1" ,"HJURP","ZNF554","DNAJA1")]
train_new_<-train_new[,-1]
head(train_new)
write.csv(train_new, file = "GSE_data_10MostImpVarDataset1.csv", row.names = FALSE)
train_new<-read.csv("GSE_data_10MostImpVarDataset1.csv",header=T)
dput(names(train_new))


head(train_new)
dim(train_new)
require(reshape)

train_new.m<-melt(train_new, id.var="Class.Label")
ggplot(train_new, aes(x=train_new$HIST1H4C,y=sqrt(train_new$HIST1H4C),fill=train_new$Class.Label))+geom_boxplot()

train_new.m
#****boxplot

ggplot(train_new.m, aes(x=variable,y=(value),fill=Class.Label))+scale_fill_brewer(palette="Dark2")+ggtitle("Top 10 Variables' Class Distribution")+geom_boxplot()

# Elaboration of heatmap with custom gradient (red - green)


ggplot(train_new.m, aes(variable, Class.Label )) +
  geom_tile(aes(fill = sqrt(value)) , color = "white") +
  scale_fill_gradient(low = "red", high = "green") +
  ylab("List of gene type ") +
  xlab("List of variables") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "value")
# Elaboration of heatmap (white - steelblue)

ggplot(train_new.m, aes(variable, Class.Label )) +
  geom_tile(aes(fill = log(value)), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of gene type ") +
  xlab("List of variables") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "value")

###*** Heatmap.2
require(gplots)
H <- train_new %>% group_by(Class.Label) %>% summarise_each(funs(mean))
H
H <- t(H)
H <- H[-1,]
H
H <- read.table("cellcycleheat.txt")
H

heatmap.2(as.matrix(H),scale="row")
heatmap.2(as.matrix(H),scale="row",col=bluered)
heatmap.2(as.matrix(H), scale="row",col=bluered,tracecol=NA)

###***
pdf("outputHetmap.pdf",width=10,height=20)
heatmap.2(as.matrix(H),scale="row",col=bluered,Colv=FALSE,tracecol=NA)###***
dev.off()

heatmap.2(as.matrix(H), tracecol=NA)

#run these if you do not have the packages installed
install.packages("ggplot2")
install.packages("Hmisc")
install.packages("dplyr")
install.packages("reshape")
install.packages("lme4")
#load libraries into R session
library(ggplot2)
library(Hmisc)
library(dplyr)
library(reshape)
library(lme4)
library(nlme)
install.packages("ggmap") # for theme_nothing
library(ggmap)
theme_set(theme_grey(base_size = 12))
#declare data and x and y aesthetics, 
#but no shapes yet
train_new_melted<-melt(train_new)

mm = melt(train_new, id='Class.Label')
mm

ggplot(mm, aes(x=Class.Label, y=value, fill=variable))+geom_bar(stat='identity', position='dodge')


mm$value<-log(mm$value)
##df <- data.frame(value = mm$value,
        ##         Group = mm$Class.Label) %>%
# factor levels need to be the opposite order of the cumulative sum of the values
##mutate(Group = factor(mm$Class.Label, levels = c("G1", "G2", "S")),
  ##     cumulative = cumsum(mm$value),
    ##   midpoint = cumulative - mm$value / 2,
      ## label = paste0(Group, " ", round(mm$value / sum(mm$value) * 100, 1), "%"))


bp<- ggplot(mm, aes(x="", y=value, fill=Class.Label))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie
# Use custom color palettes
pie + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# use brewer color palettes
pie + scale_fill_brewer(palette="Dark2")





mm = melt(train_new, id='Class.Label')
mm
p10 <- ggplot(mm, aes(x = Class.Label, y = value)) +
  geom_boxplot()
p10
p10 <- p10 + scale_x_discrete(name = "Class.Label") +
  scale_y_continuous(name = "Mean of genes coding profiles values")
p10
p10 <- p10 + scale_y_continuous(name = "Mean of genes coding profiles values",
                                breaks = seq(0, 4500, 500),
                                limits=c(0, 4500))
p10
p10 <- p10 + ggtitle("Boxplot of Cell Type Distribution ")
p10

fill <- "#4271AE"
line <- "#1F3552"

p10 <- ggplot(mm, aes(x = Class.Label, y = value)) +
  geom_boxplot(fill = fill, colour = line,
               alpha = 0.7) +
  scale_y_continuous(name = "Mean of genes coding profiles values",
                     breaks =seq(0, 4500, 500),
                     limits=c(0, 4500)) +
  scale_x_discrete(name = "Class.Label") +
  ggtitle("Boxplot of Cell Type Distribution")
p10
p10 <- p10 + theme_bw()
p10

###heat Maps
require(plyr)
dim(train_new)
df222_GSE.melted <- melt(train_new, id='Class.Label')
df222_GSE.melted <- ddply(df222_GSE.melted, .(variable), transform,rescale = log(value))
tiff("HeatMapTop10.jpeg",res=300,height=5,width=5,units="in")

(p <- ggplot(df222_GSE.melted, aes(Class.Label, variable)) + geom_tile(aes(fill = rescale),colour = "white") + scale_fill_gradient(low = "white",high = "steelblue"))
df222_GSE.s <- ddply(df222_GSE.melted, .(variable), transform,rescale = log(value))
last_plot() %+% df222_GSE.s



dev.off()

###Error barplots
# Use geom_pointrange
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))
# Use geom_line()+geom_pointrange()
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line()+
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))

###SVM ************************************************************************************
dim(df222_GSE)

###_____________________________________________



#### Model2: SVM
require("e1071")
model2 <- svm_model <- svm(Class.Label ~ ., data=train)
summary(model2)
pred2 <- predict(model2,newdata=test)
table(pred2,test.label)
confusionMatrix(data = pred2, test.label)

acc2<-(28 + 21 + 23) / nrow(test)###***0.972973
acc2
print(table(test.label))
print(table(pred2))

svm_tune <- tune(svm, train.x=x, train.y=y, kernel="linear", ranges=list(cost=10^(-2:2), gamma=2^(-2:2)))
#tune() gives us the tuned parameters, C and gamma
svm_tune
svm_model_after_tune <- svm(Class.Label ~ ., data=train,method = "C-classification", kernel="linear", cost=0.01, gamma=0.25)
require(rminer)
model<-fit(Class.Label~.,train,model="svm")
model
summary(model)

pred2.afterTune <- predict(svm_model_after_tune,test)
pred2.afterTune
system.time(predict(svm_model_after_tune,test))

table(pred2.afterTune,test.label)
confusionMatrix(data = pred2.afterTune, test.label)
acc3<-(27 + 26 + 18) / nrow(test)
acc3
table(test.label)
table(pred2.afterTune)

RocImp <- filterVarImp(x = train[, -ncol(train)], y = train$Class.Label)
RocImp
sort(RocImp[,1])
sort(RocImp[,2])
sort(RocImp[,3])

max(RocImp[,1])

install.packages("rminer")
library(rminer)
svm.imp <- Importance(model, data=train)
head(svm.imp,n=2)
sort(svm.imp$imp)

svm.imp$inputs
svm.imp$imp


VariableImportance<-Importance(model,train,method="sensv")
VariableImportance
head(VariableImportance)
L=list(runs=1,sen=t(VariableImportance$imp),sresponses=VariableImportance$sresponses)
L
mgraph(L,graph="IMP",leg=names(train),col="gray",Grid=20,cex=0.44)
###_____________________________________________________________________________________________________________________________
#### Model3:Naive Bayes:
model3 <- naiveBayes(Class.Label ~ ., data =train)
class(model3)
summary(model3)
print(model3)
pred3 <- predict(model3, newdata = test)
system.time(predict(model3,test))
table(pred3,test.label)
confusionMatrix(data = pred3, test.label)

acc4<-(28 + 20 + 20) / nrow(test)###0.9189189
acc4
print(table(test.label))
print(table(pred3))

control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"

nb_tune <- train(Class.Label~., data=train, method="nb",metric=metric, trControl=control)
nb_tune
pred3.afterTune <- predict(nb_tune, newdata = test)
system.time(predict(model3,test))
table(pred3.afterTune,test.label)
confusionMatrix(data = pred3.afterTune, test.label)
###_____________________________________________________________________________________________________________________________
#### Model4:Decision Tree:
require(rpart)
model4 <- rpart(Class.Label ~.,data=train, method="class")
class(model4)
summary(model4)
print(model4)
plot(model4)
text(model4)
pred4 <- predict(model4, test, type = "class")
system.time(predict(model4,test))
table(pred4,test.label)
confusionMatrix(data = pred4, test.label)
# prune the tree 
pfit<- prune(model4, cp=model4$cptable[which.min(model4$cptable[,"xerror"]),"CP"])
# plot the pruned tree 
plot(pfit, uniform=TRUE, 
     main="Pruned Classification Tree")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)
post(pfit, file = "ptree.ps", 
     title = "Pruned Classification Tree")
pred5 <- predict(pfit, test, type = "class")
system.time(predict(pfit,test))
table(pred5,test.label)
confusionMatrix(data = pred5, test.label)

###_____________________________________________________________________________________________________________________________
#### Model5: K- Nearest Neighbors:
train[["Class.Lable"]] = factor(train[["Class.Label"]])
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(1234)
knn_fit <- train(Class.Lable ~., data = train, method = "knn",
                 trControl=trctrl,
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
knn_fit
plot(knn_fit)
head(test)
pred6_test <- predict(knn_fit, newdata = test)
pred6_test
system.time(predict(knn_fit,test))
table(pred6_test,test.label)
confusionMatrix(data = pred6_test, test.label)###0.9459
###_____________________________________________________________________________________________________________________________
#### Model6: Neural Network:
#require(MASS)
#require(neuralnet)
#n <- names(train)
#head(train[1:20])
#m <- model.matrix( ~ Class.Label+TOP2A+UBE2C+CDK1+UBE2C.1 +HIST1H4C+CCNB1.1+CDK1.1 +CCNB1+SNRPF+ZNF554,data = train)
#head(m)
#nn <- neuralnet( Class.LabelG2 ~ TOP2A+UBE2C+CDK1+UBE2C.1 +HIST1H4C+CCNB1.1+CDK1.1 +CCNB1+SNRPF+ZNF554+HMGB2+ARRDC3+MAT2B.1+HJURP+PTBP2.1+RPL9.1+MALL.1+MRPL33.1 +DNAJA1 +SLC31A1, 
#  data=m, hidden=10, threshold=0.01
#)
# nn <- neuralnet(Class.LabelG2 ~ TOP2A+UBE2C+CDK1+UBE2C.1 +HIST1H4C+CCNB1.1+CDK1.1 +CCNB1+SNRPF+ZNF554+HMGB2+ARRDC3+MAT2B.1+HJURP+PTBP2.1+RPL9.1+MALL.1+MRPL33.1 +DNAJA1 +SLC31A1, data=m, hidden=2, err.fct="ce", linear.output=FALSE)

#plot(nn)
#head(m)
#m<-m[,-1]
#m
#compute(nn, m)


require(mlbench)
require(caret)
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3,savePred =T)
# train the LVQ model
set.seed(1234)
modelRF <- train(Class.Label~., data=reduced_Data, method="rf", trControl=control)

pred_RF <- predict(modelRF, test,type="raw")
pred_RF
system.time(predict(modelRF,test))
modelRF$pred
# train the SVM model
set.seed(1234)
modelSVM <- train(Class.Label~., data=df222_GSE, method="svmLinear", trControl=control)
pred_SVM <- predict(modelSVM, test,type="raw")
pred_SVM
system.time(predict(modelSVM,test))
confusionMatrix(data = pred_SVM, test.label)
# train the DT model
set.seed(1234)
modelDT <- train(Class.Label~., data=df222_GSE, method="rpart", trControl=control)
pred_dt <- predict(modelDT, test,type="raw")
pred_dt
system.time(predict(modelDT,test))
confusionMatrix(data = pred_dt, test.label)
# train the NB model
set.seed(1234)
modelNB <- train(Class.Label~., data=df222_GSE, method="nb", trControl=control)
pred_NB <- predict(modelNB, test,type="raw")
pred_NB
system.time(predict(modelNB,test))
confusionMatrix(data = pred_NB, test.label)
# train the KNN model
set.seed(1234)
modelKNN <- train(Class.Label~., data=df222_GSE, method="knn", trControl=control)
pred_KNN <- predict(modelKNN, test,type="raw")
pred_KNN
system.time(predict(modelANN,test))
confusionMatrix(data = pred_KNN, test.label)
# train the Model Averaged Neural Network
set.seed(1234)
modelANN <- train(Class.Label~., data=df222_GSE, method="avNNet", trControl=control)
head(test)
dim(test)
test<-test[,-114]
pred_ANN <- predict(modelANN, test,type="raw")
pred_ANN
system.time(predict(modelANN,test))
confusionMatrix(data = pred_ANN, test.label)

# collect resamples
results <- resamples(list(RF=modelRF, SVM=modelSVM, DT=modelDT,NB=modelNB,KNN=modelKNN,ANN=modelANN))
# summarize the distributions
summary(results)
Predictive_Models<-c("RF","SVM","DT","NB","KNN","ANN")
ModelAccuracy<-c(0.8929,0.9678,0.7526,0.9151,0.8744,0.8675,0.9054,0.973,0.7162,0.9459,0.9324,0.8919)
ModelError<-c(1-0.8929,1-0.9678,1-0.7526,1-0.9151,1-0.8744,1-0.8675,0,0,0,0,0,0)
dataLevels<-c("training","training","training","training","training","training","testing","testing","testing","testing","testing","testing")
"RF 0.8929  0.9054 
SVM 0.9678  0.973    
DT  0.7526  0.7162
NB  0.9151  0.9459 
KNN 0.8744  0.9324    
ANN 0.8675  0.8919"
accuracy.data <- data.frame(Predictive_Models, ModelAccuracy,dataLevels)
#accuracy.data<-cbind(accuracy.data,dataLevels)
#accuracy.data
accuracy.data
summary(accuracy.data)
str(accuracy.data)
write.csv(accuracy.data, file = "accuracy.data9.csv", row.names = FALSE)
GSE_data<-read.csv("accuracy.data9.csv",header=T)
head(GSE_data)
install.packages("Rmisc")
library(Rmisc)
dfwc <- summarySEwithin(accuracy.data, measurevar="ModelAccuracy", withinvars="dataLevels",
                        idvar="Predictive_Models", na.rm=FALSE, conf.interval=.95)
dfwc

# Make the graph with the 95% confidence interval
l<-seq(0,1,by=0.1)
pp<-ggplot(accuracy.data, aes(x= ModelAccuracy, y=ModelError, group=Predictive_Models,colour=dataLevels)) +
 
  geom_line(aes(group = Predictive_Models)) +
geom_errorbar(width=.1, aes(ymin=ModelError-ModelAccuracy, ymax=ModelAccuracy+ModelError))+
  geom_point(shape=21, size=3, fill="white")+
scale_x_continuous(limits = c(0.65, 1.80))+
coord_flip(ylim = c(-3,3), xlim = c(-2,2))
pp
########### Define the top and bottom of the errorbars
###*******************************Error Bars**********************************************************************************************

p <- ggplot(accuracy.data, aes(ModelAccuracy, Predictive_Models, colour = dataLevels))
p + geom_point() +
  geom_errorbarh(aes(xmax = ModelAccuracy + ModelError, xmin = ModelAccuracy - ModelError))
  
p + geom_point() +
  geom_errorbarh(aes(xmax =ModelAccuracy + ModelError, xmin =  ModelAccuracy - ModelError, height = .2))

###*******************************Error Bars**********************************************************************************************
####Accuracy Barplot
library(ggplot2)

library(reshape2)
library(ggplot2)
cc<-c(0,0.1,0.2,.30,.40,.50,.60,.70,.80,.90,1)
percentage<-paste(round(100*cc, 2), "%", sep="")
as.numeric(percentage)

acc_percentage<-paste(round(100*ModelAccuracy, 2), "%", sep="")
as.numeric(acc_percentage)
Accuracy_Percentage<-100*ModelAccuracy
Accuracy_Percentage
PPP<-paste(Accuracy_Percentage, "%", sep="")
PPP
Predictive_ModelsVals <- levels(accuracy.data$Predictive_Models)
accuracy.data$Predictive_Models <- as.integer(accuracy.data$Predictive_Models)
# Plot the counts
tiff("barplot2",res=300,height=5,width=5,units="in")

p<-ggplot(accuracy.data, aes(x=Predictive_Models, y=Accuracy,fill=Data)) + 
  geom_bar(stat='identity', position=position_dodge(width=1),width = 0.9) +
 # + geom_bar(width = 0.8, position = position_dodge(width = 0.9))

  scale_x_continuous(limits=c(min(accuracy.data$Predictive_Models)-1,max(accuracy.data$Predictive_Models)+1),
                     breaks=(min(accuracy.data$Predictive_Models)):(max(accuracy.data$Predictive_Models)+1), 
                     labels=c(Predictive_ModelsVals,''))+
  geom_text(aes(label=acc_percentage ), vjust=1.5, color='white', position=position_dodge(.9), size=2.2)
p  
seq<-c(0,10,20,30,40,50,60,70,80,90,100)
seq
tiff("barc",res=300,height=5,width=5,units="in")

plot <- ggplot(accuracy.data, aes(x = Predictive_Models, y = Accuracy_Percentage, fill = factor(Data))) +
  geom_bar(stat = "identity", position = "dodge", color = "grey40",width = .9) +
  scale_fill_manual(values = c("#a1d99b", "#31a354")) +
  geom_text(aes(label = acc_percentage), position = position_dodge(1),
            vjust = 1.5, color = "black", family = "Georgia",size=2.8)

plot+theme(legend.position="bottom")####*******************************************


dev.off()
getwd()
#+ scale_y_continuous(label = acc_percentage)


# Use geom_pointrange()
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")
str(training_accuracy.data)
# boxplots of results
bwplot(results,Metric="accuracy")
# dot plots of results
dotplot(results,Metric="accuracy")


results$models
results
modelANN$times
summary(modelANN)$Accuracy
#this worked 
compute(nn,as.data.frame(scaledData)[test,])

varImp(rf_gridsearch)


train_new <- train[c("Class.Label","HIST1H4C","UBE2C","UBE2C.1","CDK1" ,"CCNB1","CDK1.1","CCNB1.1" ,"TOP2A","DNAJA2","DNAJA1")]
train_new
train_new_<-train_new[,-1]
library(colorspace) # get nice colors

Class_col <- rev(rainbow_hcl(3))[as.numeric(train_new_labels)]
some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))

# scaled_iris2 <- iris2 %>% as.matrix %>% scale
# library(gplots)
d_train_new_ <- dist(train_new_) # method="man" # is a bit better
hc_train_new_ <- hclust(d_train_new_, method = "complete")
train_new_labels <- rev(levels(train[,1]))

require(dendextend)
dend <- as.dendrogram(hc_train_new_)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:150)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=3) #
require(d3heatmap)
d3heatmap::d3heatmap(as.matrix(train_new_),
                     dendrogram = "row",
                     Rowv = dend,
                     colors = "Greens",
                     # scale = "row",
                     width = 600,
                     show_grid = FALSE)
###________________________________________________________________________________________________________________________
newdata<-df222[-1]
head(newdata[1:12])
predict(resultLasso, newdata, type="response")
#resultLasso<-(glm(Class.Label ~ . ,data=df222, family=binomial))
summary(resultLasso)
head(output)
output
tail(output)
table(output[4]<0.05)
dim(output)
names(resultLasso)
resultLasso$coefficients
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
### Junk
#head(GSE_data_2[1:4],n=45)
tail(GSE_data[1:12],n=5)
dim(GSE_data_new)
GSE_data_2
head(GSE_data_2)
i
omitted_col
#GSE_data2<-NULL
as.data.frame(GSE_data2)
str(GSE_data2)
dim(GSE_data_2)
GSE_data2<-NULL
head(GSE_data2[1:12],n=12)
tail(GSE_data2[1:12],n=12)
as.data.frame(GSE_data2)
G1<-subset(GSE_data,select=ATP13A2)
G1
cbind(GSE_data2[2],G1)
GSE_data2$ATP13A2<-NULL
GSE_data2<-cbind(GSE_data2[2],G1)

S<-nrow(subset(GSE_data[92:171,], GSE_data[92:171,][,2]!="0"))
S
G2<-nrow(subset(GSE_data[172:274,], GSE_data[172:274,][,2]!="0"))
G2
G1<-nrow(subset(GSE_data[1:92,],GSE_data[1:92,][,2]!="0"))
G1
S<-nrow(subset(GSE_data[93:172,], GSE_data[93:172,][,2]!="0"))
S
G2<-nrow(subset(GSE_data[173:248,], GSE_data[173:248,][,2]!="0"))
G2

#for(i in 1:130){
 # print(i)
  #cola <- paste('col', i, sep= '')
  #df[[cola]] <- ifelse(x == i, 1, 0)
#}
