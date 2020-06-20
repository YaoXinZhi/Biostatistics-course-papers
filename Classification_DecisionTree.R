# e1071 -- svm
# randomForest -- RandomForest
# glmnet -- lasso ridge


#pkgs <- c('rpart', 'rpart.plot', 'party', 'randomForest', 'e1071')
#install.packages(pkgs, depend=TRUE)
install.packages('rpart')


setwd('Desktop/Biostatistics')
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
  
# DecisionTree
library('rpart')
dtree <- rpart(class~., data=df.train, method='class',
               parms=list(split="information"))
dtree$cptable

print('summary: ')
summary(dtree)

plotcp(dtree)
# pruned
dtree.pruned <- prune(dtree, cp=0.01253)
# plot
library('rpart.plot')
?prp
prp(dtree.pruned, type=2, extra=104,
    fallen.leaves=TRUE, main='Decision Tree')
# prediction
dtree.pred <- predict(dtree.pruned, df.validate, type='class')
dtree.perf <- table(df.validate$class, dtree.pred,
                    dnn=c('Actual', 'Predicted'))
dtree.perf

performance(dtree.perf)

