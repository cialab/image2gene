### Examine the accuracy of utilizing the top genes (all 15) ###

set.seed(1)
require(xgboost)
require(caret)
require(parallel)
require(parallelMap)
require(data.table)
parallelStartSocket(cpus = detectCores(), show.info = FALSE)
require(mlr)
configureMlr(show.info = FALSE)
require(e1071)

# Top 15 genes (in reverse order)
# 22922,5814,6598,201,6651,6649,16896,12114,6597,7001,4337,18325,8187,8728,8185
topk=c("x74145","x16178","x17395","x100038909","x17750","x17748","x546644","x245195","x17394","x18407","x12765","x64381","x20717","x213002","x20716")
topk=c("x74145","x16178","x17395","x100038909","x17750","x17748","x546644","x245195","x17394","x18407","x12765","x64381","x20717","x213002")
topk=rev(topk)

# Read data
data = read.table('../data.csv', header = T, sep = ',')


for (asdf in seq(1,100)){
for (k in seq(1,15)){
  # Store stuff here
  preds = rep(0,90)
  gts = rep(0,90)
  # LOO
  for (i in seq(1,90)){
    
    # I love this minus syntax
    train <- data[-i, c("Class",topk[1:k])]
    test <- data[i,c("Class",topk[1:k])]
    
    # Change to two class
    train$Class=as.factor(as.numeric(train$Class=="SS"))
    test$Class=factor(as.numeric(test$Class=="SS"), levels = c(0,1))
    traintask <- makeClassifTask(data = train, target = "Class", positive = 1)
    testtask <- makeClassifTask(data = test, target = "Class", check.data = FALSE, fixup.data = "no", positive = 1)
    
    # Learner
    lrn <- makeLearner("classif.xgboost",predict.type = "response")
    lrn$par.vals <- list(
      objective="binary:logistic",
      eval_metric="logloss"
    )
    
    # Parameter space
    params <- makeParamSet(
      makeDiscreteParam("booster",values = c("gbtree")),
      makeIntegerParam("nrounds", lower = 1, upper = 100),
      makeIntegerParam("max_depth",lower = 3, upper = 10),
      makeNumericParam("eta", lower = 0.01, upper = 0.5),
      makeNumericParam("min_child_weight",lower = 1,upper = 10),
      makeNumericParam("subsample",lower = 0.5,upper = 1),
      makeNumericParam("colsample_bytree",lower = 0.5,upper = 1),
      makeNumericParam("lambda", lower = -1, upper = 0, trafo = function(x) 10^x))
    
    # Resampling strategy
    rdesc <- makeResampleDesc("CV",stratify = T, iters = 10)
    
    # Search strategy
    ctrl <- makeTuneControlRandom(maxit = 100)
    
    # Parameter tuning
    mytune <- tuneParams(
      learner = lrn,
      task = traintask,
      resampling = rdesc,
      measures = acc,
      par.set = params,
      control = ctrl,
      show.info = FALSE)
    
    # Set hyperparameters
    lrn_tune <- setHyperPars(lrn, par.vals = mytune$x)
    
    # Train optimal model
    xgmodel <- train(learner = lrn_tune, task = traintask)
    
    # Predict training set
    #pred <- predict(xgmodel,traintask)
    #a = confusionMatrix(pred$data$response,pred$data$truth, dnn = c("Pr","GT"), positive = "1")
    
    # Predict testing set
    pred <- predict(xgmodel,testtask)
    preds[i]=as.numeric(as.character(pred$data$response))
    gts[i]=as.numeric(as.character(pred$data$truth))
    
    # Train on whole training set
    #params = list(booster="gbtree",
    #              max_depth=mytune$x$max_depth,
    #              eta=mytune$x$eta,
    #              min_child_weight=mytune$x$min_child_weight,
    #              subsample=mytune$x$subsample,
    #              colsampe_bytree=mytune$x$colsample_bytree,
    #              lambda=mytune$x$lambda,
    #              objective="binary:logistic")
    #dtrain = xgb.DMatrix(label=as.numeric(train$Class)-1, data = as.matrix(train[,-1]))
    #mdl=xgb.train(params, dtrain, nrounds = mytune$x$nrounds,  eval_metric = "logloss")
    
    # Importance
    #imp<-xgb.importance(mdl$feature_names,model = mdl)
    #write.table(imp, paste('fourthTry/features',toString(i),'.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
    #write.table(mytune$x, paste('fourthTry/params',toString(i),'.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
  }
  a=confusionMatrix(as.factor(preds),as.factor(gts), dnn = c("Pr","GT"), positive = "1")
  write.table(a$table, paste('topk_minus1/te',k,'_',asdf,'.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
}
}
parallelStop()