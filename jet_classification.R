library(ggplot2)
library(plyr)
library(dplyr)
library(caret)
library(cvTools)
library(ROCR)
library(class)
library(nnet)
# loading in my jet data
# the bulk of my data cleaning is done before this R script

jet_data <- read.csv("C:/Users/Student/Desktop/DS 4559 2016/final/jet_data.txt", sep="")

# jet_data$source represent the particle that caused the jet
# each number represents a different type of particle
# negative numbers represent anti-particles
# I will now group them into
# 1 = light quarks (up, down, strange)
# 2 = top
# 3 = gluons 
jet_data$source <- abs(jet_data$source)
jet_data <- jet_data[!jet_data$source == 4,]
jet_data$source[jet_data$source < 3.5] <- 1
jet_data$source[jet_data$source == 6] <- 2
jet_data$source[jet_data$source == 21] <- 3

#cleaning measure to filter out bad jets
jd<- jet_data[jet_data$numCS > 10 & (jet_data$pt > 1000 & jet_data$pt < 300000),]



#looking at pt distributions
jd_1 <- data.frame(jd[jd$source == 1,])
jd_2 <- data.frame(jd[jd$source == 2,])
jd_3 <- data.frame(jd[jd$source == 3,])

hist(jd_3$pt, border = "red", freq = FALSE, breaks = 50, xlim = c(0,300000),main = "Transverse Momentum", xlab = "pT")
hist(jd_2$pt, add=T, border = "blue", freq = FALSE, breaks = 100, xlim = c(0,300000))
hist(jd_1$pt, add=T, border = "green", freq = FALSE, breaks = 50, xlim = c(0,300000))
legend(180000, .000008, c("gluons","top quarks", "light quarks"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("red","blue","green"))

hist(jd_3$m, border = "red", freq = FALSE, breaks = 50, xlim = c(0,150000),main = "Invarient Mass", xlab = "mass")
hist(jd_2$m, add=T, border = "blue", freq = FALSE, breaks = 100, xlim = c(0,150000))
hist(jd_1$m, add=T, border = "green", freq = FALSE, breaks = 50, xlim = c(0,150000))
legend(90000, .00003, c("gluons","light quarks", "top quarks"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("red","blue","green"))

hist(jd_1$numCS, border = "red", freq = FALSE, breaks = 30, xlim = c(10,120),main = "NUmber of Constituents", xlab = "Rapidity Difference")
hist(jd_2$numCS, add=T, border = "blue", freq = FALSE, breaks = 30, xlim = c(10,120))
hist(jd_3$numCS, add=T, border = "green", freq = FALSE, breaks = 30, xlim = c(10,120))
legend(75, .03, c("light","top quarks", "gluons"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("red","blue","green"))


jd <- subset(jd, select = c(source, pt, mt, m, numCS, E, avgR2Weighted, avgRapDiff, avgR2))
jd_rand <- jd[order(runif(nrow(jd))),]
jd_rand$source <- as.factor(jd_rand$source)
normalize <- function(x) {
  if(is.numeric(x)){
    return((x-min(x))/(max(x) - min(x)))
  }
  return(x)
}
jd_norm <- as.data.frame(lapply(jd_rand, normalize))

library(randomForest)
library(C50)
library(RWeka)
set.seed(1)
#variable_importance <- jd_norm
#variable_importance$random <- runif(nrow(jd_norm))
#jd_m <- randomForest(source ~ ., variable_importance)
#varImpPlot(jd_m)


## C5.0/JRip
## Divide data into j folds
j = 10
folds <- cvFolds(NROW(jd_rand), K=j)
dataset <- jd_norm
dataset$holdoutpred <- rep(0,nrow(dataset))
dataset$holdoutprob <- rep(0,nrow(dataset))

#code is in here for JRip, knn, or C5.0. its up to the user to decide which to use
for(i in 1:j){
  train_data <- dataset[folds$subsets[folds$which != i], ] #Set the training set
  ideal <- class.ind(train_data$source)
  validation_data <- dataset[folds$subsets[folds$which == i], ] #Set the validation set
  #newmod <- nnet(train_data[,-1],ideal, size=10,softmax = TRUE)
  #newmod <- JRip(source ~.,data=train_data)
  #newpred <- knn(train=train_data,test=validation_data,cl=train_data[,1],k=10)
  newmod <- C5.0(source ~ .,data=train_data) #Get your new linear model (just fit on the train data)
  
  newpred <- predict(newmod,newdata=validation_data[,-1]) #Get the predicitons for the validation set (from the model just fit on the train data)
  #newprob <- predict(newmod, newdata=validation_data[,-1],type="prob")[,2]
  dataset[folds$subsets[folds$which == i], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use.
  #dataset[folds$subsets[folds$which == i], ]$holdoutprob <- newprob
}

pred <- as.factor(dataset$holdoutpred)
cm<-confusionMatrix(pred, jd_test$source)
cm_table <- unname(cm$table)
colnames(cm_table) = c("light","top", "gluon")
rownames(cm_table) = c("light","top", "gluon")
for(i in 1:3){
  cm_table[,i] <- cm_table[,i]/sum(cm_table[,i])
}
#plotting the confusion matrix
confusion <- as.data.frame(t(as.table(data.matrix(cm_table))))

plot <- ggplot(confusion)
plot  +geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="Actual Class") + scale_y_discrete(name="Predicted Class") + scale_fill_gradient(breaks=seq(from=-.5, to=4, by=.2)) +  labs(fill="Normalized\nFrequency") + geom_text(aes(x = Var1, y = Var2, label = sprintf("%3.2f", Freq)), vjust = 1)

