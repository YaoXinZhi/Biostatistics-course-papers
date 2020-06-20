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


## RandomForest
library(randomForest)
# build forest
print('running')
fit.forest <- randomForest(class~., data=df.train,
                           na.action=na.roughfix,
                           importance=TRUE)
fit.forest
# give variable importance
#importance(fit.forest, type=2)
# classify sample points outside the training data
forest.pred <- predict(fit.forest, df.validate)
length(forest.pred)
length(df.validate$class)
forest.perf <- table(df.validate$class, forest.pred,
                     dnn=c('Actual', 'predicted'))
forest.perf


performance(forest.perf)

