set.seed(1)

require(xgboost)
require(caret)
require(parallel)
require(parallelMap)
require(data.table)
parallelStartSocket(cpus = detectCores(), show.info = FALSE)
require(mlr)
configureMlr(show.info = FALSE)

# Read data
data0 = read.table('../data.csv', header = T, sep = ',')
folds = read.table('folds.txt', header = F, sep = ',', fill = TRUE)
ss=c()
for (i in seq(1,10)){ ss = c(ss,folds[i,])}
ss=as.numeric(ss)
ss=ss[!is.na(ss)]
data=data0[data0$MouseNum %in% ss,]
data33 = read.table('../data33.csv', header = T, sep = ',')
wrun="run_handpicked_test"

# 8185,8728,8187,18325,4337
genes = c("x20716","x213002","x20717","x64381","x12765")

# Training the model
train <- data[,c("Class",genes)]
train$Class=as.factor(as.numeric(train$Class=="SS"))
traintask <- makeClassifTask(data = train, target = "Class", positive = 1)
lrn <- makeLearner("classif.xgboost",predict.type = "response")
lrn$par.vals <- list(
  objective="binary:logistic",
  eval_metric="logloss")
params <- makeParamSet(
  makeDiscreteParam("booster",values = c("gbtree")),
  makeIntegerParam("nrounds", lower = 1, upper = 100),
  makeIntegerParam("max_depth",lower = 3, upper = 10),
  makeNumericParam("eta", lower = 0.01, upper = 0.5),
  makeNumericParam("min_child_weight",lower = 1,upper = 10),
  makeNumericParam("subsample",lower = 0.5,upper = 1),
  makeNumericParam("colsample_bytree",lower = 0.5,upper = 1),
  makeNumericParam("lambda", lower = -1, upper = 0, trafo = function(x) 10^x))
rdesc <- makeResampleDesc("CV",stratify = T, iters = 10)
ctrl <- makeTuneControlRandom(maxit = 100)
mytune <- tuneParams(
  learner = lrn,
  task = traintask,
  resampling = rdesc,
  measures = acc,
  par.set = params,
  control = ctrl,
  show.info = FALSE)
lrn_tune <- setHyperPars(lrn, par.vals = mytune$x)
params = list(booster="gbtree",
              max_depth=mytune$x$max_depth,
              eta=mytune$x$eta,
              min_child_weight=mytune$x$min_child_weight,
              subsample=mytune$x$subsample,
              colsampe_bytree=mytune$x$colsample_bytree,
              lambda=mytune$x$lambda,
              objective="binary:logistic")
dtrain = xgb.DMatrix(label=as.numeric(train$Class)-1, data = as.matrix(train[,-1]))
mdl=xgb.train(params, dtrain, nrounds = mytune$x$nrounds,  eval_metric = "logloss")

# Folds
for (f in seq(1,10)){
  # Nums
  for (num_instances in c(2500,10000)){#00000)){
    # Sampling
    for (sampling in c("random","uniform","segments")){
      if (sampling=="segments"){
        prefix=paste(sampling,'_',as.integer(num_instances),'_2',sep='')
      }
      else{
        prefix=paste(sampling,'_',as.integer(num_instances),sep='')
      }
      
      # Predicted gene values
      preds=read.table(text = gsub("\t"," ", 
                                   gsub("\\]"," ", 
                                        gsub('\\['," ",readLines(paste('X:/python/AttentionDeepMIL/image2gene/revision/run_handpicked_test/',
                                                                       prefix,'_',f-1,'_test.txt', sep = ''))))))
  
      # Get their classes :(
      asdf=numeric(dim(preds)[1])
      for (i in seq(1,dim(preds)[1])){
        if (data33$Class[data33$MouseNum==preds[i,1]]=="SS")
          asdf[i]=1
      }
  
      # Ground truth
      test2gt <- preds[,c(2,3,4,5,6)]
      colnames(test2gt) <- genes
      test2gtc <- asdf
      dtest2gt = xgb.DMatrix(label=test2gtc, data = as.matrix(test2gt))
      pr=round(predict(mdl,dtest2gt))
      write.table(t(pr),paste('on_predictions_test/',prefix,'_gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
      
      # Predicted values
      test2pr <- preds[,c(7,8,9,10,11)]
      test2n <- preds[,c(1)]
      colnames(test2pr) <- genes
      test2prc <- asdf
      dtest2pr = xgb.DMatrix(label=test2prc, data = as.matrix(test2pr))
      pr=round(predict(mdl,dtest2pr))
      write.table(t(pr),paste('on_predictions_test/',prefix,'_pr.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
    }
  }
}
write.table(t(preds$V1),paste('on_predictions_test/gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
write.table(t(asdf),paste('on_predictions_test/gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)

parallelStop()