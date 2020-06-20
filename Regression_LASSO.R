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
sgRNA.1130<-read.csv(file="data/Regression/sgRNA.1130.csv",header=T)
sgRNA.1130.names<-as.character(sgRNA.1130[,1])
sgRNA.top20.id<-sgRNA.1130.names[1:20]
rownames(x)<-sgRNA.1130.names
y<-sgRNA.1130[,5]
names(y)<-sgRNA.1130.names

#length(X)
#length(sgRNA.1130.names)
# Construct the training groups and testing groups
train.groups.id<-list()
test.groups.id<-list()
x.train<-list()
x.test<-list()
y.train<-list()
y.test<-list()
y.test.total<-c()
tmp<-c()
for(i in 1:10){
  test.groups.id[[i]]<-sgRNA.1130.names[which(sgRNA.1130[,6]==i)]
  train.groups.id[[i]]<-sgRNA.1130.names[which(sgRNA.1130[,6]!=i)]
  x.train[[i]]<-x[train.groups.id[[i]],]
  x.test[[i]]<-x[test.groups.id[[i]],]
  y.train[[i]]<-y[train.groups.id[[i]]]
  y.test[[i]]<-y[test.groups.id[[i]]]
  tmp<-c(tmp,test.groups.id[[i]])
  y.test.total<-c(y.test.total,y.test[[i]])
}
names(y.test.total)<-tmp

length(sgRNA.1130)
##Perform Lasso
library(glmnet)
sum_y_top20<-rep(0,20)
sum_y_total<-rep(0,1130)
for(m in 1:10){
  cat("m=",m,'\t')
  y.test.pred.total<-c()
  for(i in 1:10){
    Lasso.cv<-cv.glmnet(x.train[[i]],y.train[[i]])
    bestlambda<-Lasso.cv$lambda.min
    Lasso.model<-Lasso.cv$glmnet.fit
    y.test.pred.Lasso<-predict(Lasso.model,newx=x.test[[i]],s=bestlambda)
    y.test.pred.total<-c(y.test.pred.total,y.test.pred.Lasso)
  }
  names(y.test.pred.total)<-tmp
  sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
  sum_y_total <- y.test.pred.total+sum_y_total
}

y.test.pred.top20 <- sum_y_top20/10
y.test.pred.total.Lasso.mean <- sum_y_total/10
PCC.Lasso.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Lasso.7mer.all <- cor(y.test.total,y.test.pred.total.Lasso.mean)
save(y.test.total,y.test.pred.total.Lasso.mean,PCC.Lasso.7mer.top20,PCC.Lasso.7mer.all,file="PCC.Lasso.200.6mer.Rdata")
# fixme:
print('summary')
summary(Lasso.cv)
print('residuals')
residuals(Lasso.cv)


