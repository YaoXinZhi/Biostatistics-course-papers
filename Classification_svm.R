  # e1071 -- svm
  # randomForest -- RandomForest
  # glmnet -- lasso ridge
  
  
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


## SVM


#print('tune svm')

#tuned <- tune.svm(class~., data=df.train,
#                  gamma=10^(-6:1),
#                  cost=10^(-10:10))

#tuned

# cost 1, gamma 0.0002441406
library('e1071')
fit.svm <- svm(class~., data=df.train,
               gamma=0.0002441406, cost=1)
fit.svm
svm.pred <- predict(fit.svm, na.omit(df.validate))
svm.perf <- table(na.omit(df.validate)$class,
                  svm.pred, dnn=c('Actual', 'Predicted'))
svm.perf

#performance(svm.perf)
