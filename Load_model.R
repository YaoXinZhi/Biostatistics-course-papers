# Read Lasso model

library(glmnet)
Lasso.model <- load('result/PCC.Lasso.50.6mer.Rdata')
Lasso.model <- load('result/PCC.Lasso.100.6mer.Rdata')
Lasso.model <- load('result/PCC.Lasso.200.6mer.Rdata')

PCC.Lasso.7mer.all
PCC.Lasso.7mer.top20
y.test.pred.total.Lasso.mean
y.test.total

write.table(y.test.pred.total.Lasso.mean,file='result/lasso_pred/lasso.200.pred.txt',sep ="\t")
write.table(y.test.total, file='result/lasso_pred/lasso.200.true.txt',sep='\t')

# Read Redige model
ridge.model <- load('result/PCC.Ridge.50.6mer.Rdata')
Lasso.model <- load('result/PCC.Ridge.100.6mer.Rdata')
Lasso.model <- load('result/PCC.Ridge.200.6mer.Rdata')

PCC.Ridge.7mer.all
PCC.Ridge.7mer.top20
y.test.pred.total.Ridge.mean
y.test.total
  
write.table(y.test.pred.total.Ridge.mean,file='result/ridge_pred//ridge.200.pred.txt',sep ="\t")
write.table(y.test.total, file='result/ridge_pred//ridge.200.true.txt',sep='\t')




