library(randomForest)
library(datasets)
library(caret)
library(ROCR)
library(pROC)
library(ggRandomForests)



clus2<-read.csv("all_tissues_Ruth_2022.csv", sep = ",",header = T, row.names = 1)
clus2$Tissue<-as.factor(clus2$Tissue)
clus2$Treatment<-as.factor(clus2$Treatment)
clus2$Location<-as.factor(clus2$Location)


clus2[is.na(clus2)] <- 0 ##Better - replaces NA with 0
clus2<-clus2[c("Tissue","AQP4","ATP4A","CCKBR","CHIA","CLIC6","GAST","GHRL","GIF","HDC","HRH2","IL8","KCNQ1","MBOAT","MUC1","MUC2","MUC5AC","MUC6","OLFM4","PCSK1","PGA5","PIGR","SST")]

#another way#
clus5 <- clus2[clus2$Tissue %in% c("Pyloric", "Cardiac"),]# Multiple factor levels
clus5 <- droplevels(clus5)  #takes out extra factor levels
summary(clus5)

clus3 <- clus2[clus2$Tissue %in% c("Fundic", "Cardiac"),]# Multiple factor levels
clus3 <- droplevels(clus3)  #takes out extra factor levels
summary(clus3)

clus4 <- clus2[clus2$Tissue %in% c("Fundic", "Pyloric"),]# Multiple factor levels
clus4 <- droplevels(clus4)  #takes out extra factor levels
summary(clus4)


summary(clus2)
#Select test and train group#
set.seed(222)
ind <- sample(2, nrow(clus2), replace = TRUE, prob = c(0.6, 0.4))
train <- clus2[ind==1,]
test <- clus2[ind==2,]
#random forest test group#
rf <- randomForest(Tissue~., data=train, proximity=TRUE, importance=TRUE)

###Plot the error rate for each tissue with legend##
plot(rf, legend("topright", legend=unique(train$Tissue), col=unique(as.numeric(train$Tissue)), pch=19))


varImpPlot(rf) 
#confusion matrix training#
p1 <- predict(rf, train)
confusionMatrix(p1, train$Tissue)
#confusion matrix test#
p2 <- predict(rf, test)
confusionMatrix(p2, test$Tissue)
as.table(p2)
##See how the predictions have been assigned to each sample##
write.csv(p2,file = "predictions_test.csv")

plot(rf)
MDSplot(rf, train$Tissue)
summary(test)
summary(train)

write.csv(confusion$mode, file="confusion.test.csv")
          
t <- tuneRF(train[,-5], train[,5],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)
t
plot(t)

hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance#
varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
TOP_10_RF<-clus2[c("Tissue","MUC5AC","HDC","GIF","HRH2","MUC2","GAST","ATP4A","CHIA","CCKBR","PGA5")]

importance(rf, cond)
imp_caret<-varImp(rf,conditional=TRUE)
varim

confusionMatrix(reference = test$Tissue, data = predicted, mode='everything', positive='MM')

#### cardiac pyloric###
#Select test and train group#
set.seed(222)
ind4 <- sample(2, nrow(clus4), replace = TRUE, prob = c(0.6, 0.4))
train4 <- clus4[ind4==1,]
test4 <- clus4[ind4==2,]
#random forest test group#
rf4 <- randomForest(Tissue~., data=train4, proximity=TRUE, importance=TRUE)
plot(rf4)
varImpPlot(rf4,main = "Pyloric v Fundic- Variable Importance") 
rf4

p1_4 <- predict(rf4, train4)
confusionMatrix(p1_4, train4$Tissue)
#confusion matrix test#
p2_4 <- predict(rf4, test4)
confusion_Test4<-confusionMatrix(p2_4, test4$Tissue)
confusion_Test4
varImpPlot(rf4)

#plot random forest elements#
ggrf2<-rfsrc(Tissue~., data = clus2, importance = TRUE)

summary(ggrf2)
ggrf2$event.info
ggrf2
##ggrandomforest##
# Plot the VIMP rankings of independent variables.
plot(gg_vimp(ggrf2))


gg_md <- gg_minimal_depth(ggrf2)
# plot the object
plot(gg_md)

gg_md_vimp<-gg_minimal_vimp(ggrf2)
plot(gg_md_vimp)



dev.off
#ROC curve#
pred1=predict(rf5,type = "prob")
library(ROCR)
perf = prediction(pred1[,1], test5$Tissue)
# 1. Area under curve
auc = performance(perf, "auc")
auc
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
plot(pred3,main="ROC Curve for Cardiac v Pyloric",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
