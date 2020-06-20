
# First you need to convert the loaction file to fasta file
# Then convert the fasta file to the feature file through the kmer.py file
# Finally, calculate the Z score through RF.model.Rdata

# e1071 -- svm
# randomForest -- RandomForest
# glmnet -- lasso ridge

# library('randomForest')
# x = as.matrix(read.table('example_7mer.txt'))
# load('RF.model.Rdata')
# y_pred = predict(RF.model, x)

#pkgs <- c('glmnet', 'rpart.plot', 'kernlab', 'randomForest', 'e1071')
#install.packages(pkgs, depend=TRUE)

#setwd('Desktop/Biostatistics')
# Load the data of input and output
x<-as.matrix(read.table(file="data/Regression/seq/p53.1130sgRNA_200bp_6mer.txt"))
sgRNA<-read.csv(file="data/Regression/sgRNA.1130.csv",header=T)
sgRNA.names<-as.character(sgRNA[,1])
sgRNA.top20.id<-sgRNA.names[1:20]
rownames(x)<-sgRNA.names
y<-sgRNA[,5]
names(y)<-sgRNA.names

#length(X)
#length(sgRNA.names)
# Construct the training groups and testing groups
train.groups.id<-list()
test.groups.id<-list()
x.train<-list()
x.test<-list()
y.train<-list()
y.test<-list()
y.test.total<-c()
temp<-c()
for(i in 1:10){
  test.groups.id[[i]]<-sgRNA.names[which(sgRNA[,6]==i)]
  train.groups.id[[i]]<-sgRNA.names[which(sgRNA[,6]!=i)]
  x.train[[i]]<-x[train.groups.id[[i]],]
  x.test[[i]]<-x[test.groups.id[[i]],]
  y.train[[i]]<-y[train.groups.id[[i]]]
  y.test[[i]]<-y[test.groups.id[[i]]]
  temp<-c(temp,test.groups.id[[i]])
  y.test.total<-c(y.test.total,y.test[[i]])
}
names(y.test.total)<-temp

length(sgRNA)

# #setwd('Desktop/Biostatistics')
#   set.seed(1234)
#   # read data
#   data <- as.data.frame(read.table('data/Regression/regression.feature.txt', sep=',', header=TRUE))  
#   length(data)
#   print('dim(data): ')
#   dim(data)
#   train <- sample(nrow(data), 0.7*nrow(data))
#   df.train <- data[train,]
#   df.validate <- data[-train,]
#   dim(df.train)
#   dim(df.validate)

fit.top<-rep(0,20)
sum_y_total<-rep(0,1130)
for(m in 1:10){
  cat("m=",m,'\t')
  fit.pred.all<-c()
  for(i in 1:10){
    linear.fit <- lm(x.train[[i]],y.train[[i]])
    bestlambda<-linear.fit$lambda.min
    linear.model<-linear.fit$glmnet.fit
    fit.pred<-predict(linear.model,newx=x.test[[i]],s=bestlambda)
    fit.pred.all<-c(fit.pred.all,fit.pred)
  }
  names(fit.pred.all)<-temp
  fit.top <- fit.pred.all[sgRNA.top20.id] + fit.top
  sum_y_total <- fit.pred.all+sum_y_total
}

y.test.pred.top20 <- fit.top/10
fit.pred.all.Linear.mean <- sum_y_total/10
PCC.Linear.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Linear.7mer.all <- cor(y.test.total,fit.pred.all.Linear.mean)
save(y.test.total,fit.pred.all.Linear.mean,PCC.Linear.7mer.top20,PCC.Linear.7mer.all,file="PCC.linear.50.6kmer.Rdata")
# fixme:
print('summary')
summary(linear.fit)
print('residuals')
residuals(linear.fit)

  


