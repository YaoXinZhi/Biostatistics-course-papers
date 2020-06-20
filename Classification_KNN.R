# e1071 -- svm
# randomForest -- RandomForest
# glmnet -- lasso ridge

#install.packages('kknn')



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

library('kknn')
#?cv.kknn
#print('tune knn')
# 10-flod cross validation
#gold.knn <- cv.kknn(class~., data=df.train, kcv=10)
#glod.knn <- train.kknn(class~., data=df.train, kernel=c("rectangular", "triangular", "epanechnikov", "optimal"),distance=2,scale=T)
# error rate
#gold.knn$MISCLASS
#print('gold.knn')
#gold.knn
# best k ==2

fit.knn <- kknn(class~., df.train, df.validate, k=2, scale=T,distance=1, kernel='rectangular')
summary(fit.knn)
fit <- fitted(fit.knn)

table(fit, df.validate$class)

# knn.pred <- predict(fit.knn, na.omit(df.validate))
# knn.perf <- table(na.omit(df.validate)$class,
#                   knn.pred, dnn=c('Actual', 'Predicted'))
# knn.perf



