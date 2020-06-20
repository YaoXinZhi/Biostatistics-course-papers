# e1071 -- svm
# randomForest -- RandomForest
# glmnet -- lasso ridge


#pkgs <- c('rpart', 'rpart.plot', 'party', 'randomForest', 'e1071')
#install.packages(pkgs, depend=TRUE)

#setwd('Desktop/Biostatistics')
# read data
data <- as.data.frame(read.table('data/Classification/ABI5-total.txt', sep=',', header=TRUE))  
  data$class <- factor(data$class, levels=c(0,1),
                       labels=c('pos', 'neg'))
table(data$class)
set.seed(1234)
train <- sample(nrow(data), 0.7*nrow(data))
df.train <- data[train,]
df.validate <- data[-train,]
table(df.train$class)
table(df.validate$class)
  

# Logstic Regression
fit.logit <- glm(class~., data=df.train, family=binomial())
summary(fit.logit)
prob <- predict(fit.logit, df.validate, type='response')
logit.pred <- factor(prob > .5, levels=c(FALSE, TRUE),
                     labels=c('pos', 'neg'))
logit.perf <- table(df.validate$class, logit.pred,
                    dnn=c('Actual', 'Predicted'))
logit.perf

