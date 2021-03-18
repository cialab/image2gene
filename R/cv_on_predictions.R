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
wrun="run_handpicked"

# 8185,8728,8187,18325,4337
genes = c("x20716","x213002","x20717","x64381","x12765")


# Folds
gt=rep(0,77)
te=rep(0,77)
teSanity=rep(0,77)
ns=rep(0,77)
sss=1
for (f in seq(1,10)){
  ss=as.numeric(folds[f,])
  ss=ss[!is.na(ss)]
  for (s in seq(1,length(ss))){
    ns[sss]=folds[f,s]
    # What even is this language
    #train <- data[!data$MouseNum %in% c(unlist(folds[f,],use.names = FALSE)), c("Class",genes)]
    #test <- data[data$MouseNum %in% c(unlist(folds[f,],use.names = FALSE)), c("Class",genes)]
    train <- data[!data$MouseNum %in% c(folds[f,s]), c("Class",genes)]
    test <- data[data$MouseNum %in% c(folds[f,s]), c("Class",genes)]
    
    # Optimization stuff
    train$Class=as.factor(as.numeric(train$Class=="SS"))
    test$Class=factor(as.numeric(test$Class=="SS"), levels = c(0,1))
    gt[sss]=as.numeric(as.character(test$Class[1]))
    
    traintask <- makeClassifTask(data = train, target = "Class", positive = 1)
    testtask <- makeClassifTask(data = test, target = "Class", check.data = FALSE, fixup.data = "no", positive = 1)
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
    
    # Retrain with optimal parameters training set
    params = list(booster="gbtree",
                  max_depth=mytune$x$max_depth,
                  eta=mytune$x$eta,
                  min_child_weight=mytune$x$min_child_weight,
                  subsample=mytune$x$subsample,
                  colsampe_bytree=mytune$x$colsample_bytree,
                  lambda=mytune$x$lambda,
                  objective="binary:logistic")
    dtrain = xgb.DMatrix(label=as.numeric(train$Class)-1, data = as.matrix(train[,-1]))
    dtest = xgb.DMatrix(label=as.numeric(test$Class)-1, data = as.matrix(test[,-1]))
    mdl=xgb.train(params, dtrain, nrounds = mytune$x$nrounds,  eval_metric = "logloss")
    
    # Eval on ground truth gene values
    #pr=round(predict(mdl,dtrain))
    #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(train$Class, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
    #write.table(aa$table, paste('on_predictions/',prefix,'_','tr.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
    pr=round(predict(mdl,dtest))
    te[sss]=pr
    #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(test$Class, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
    #write.table(aa$table, paste('on_predictions/',prefix,'_','te.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
        
    # Nums
    for (num_instances in c(2500)){#,10000,100000)){
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
                                          gsub('\\['," ",readLines(paste('X:/python/AttentionDeepMIL/image2gene/revision/run_handpicked/',
                                                                         prefix,'_',f-1,'_test.txt', sep = ''))))))
        
        # Get their classes :(
        asdf=numeric(dim(preds)[1])
        for (i in seq(1,dim(preds)[1])){
          if (data$Class[data$MouseNum==preds[i,1]]=="SS")
            asdf[i]=1
        }
        
        # Sanity check
        #train2gt <- preds[!preds[,1] %in% c(unlist(folds[f,],use.names = FALSE)), c(2,3,4,5,6)]
        #test2gt <- preds[preds[,1] %in% c(unlist(folds[f,],use.names = FALSE)),c(2,3,4,5,6)]
        test2gt <- preds[preds[,1] %in% c(folds[f,s]),c(2,3,4,5,6)]
        #colnames(train2gt) <- genes
        colnames(test2gt) <- genes
        #train2gtc <- asdf[!preds[,1] %in% c(unlist(folds[f,],use.names = FALSE))]
        #test2gtc <- asdf[preds[,1] %in% c(unlist(folds[f,],use.names = FALSE))]
        test2gtc <- asdf[preds[,1] %in% c(folds[f,s])]
        #dtrain2gt = xgb.DMatrix(label=train2gtc, data = as.matrix(train2gt))
        dtest2gt = xgb.DMatrix(label=test2gtc, data = as.matrix(test2gt))
        # pr=round(predict(mdl,dtrain2gt))
        #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(train2gtc, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
        #write.table(aa$table, paste('mil_cv/',wrun,'/trSanity.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
        #pr=round(predict(mdl,dtest2gt))
        #teSanity[sss]=pr
        #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(test2gtc, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
        #write.table(aa$table, paste('on_predictions/',prefix,'_','teSanity.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
        
        # Predicted values
        #train2pr <- preds[!preds[,1] %in% c(unlist(folds[f,],use.names = FALSE)), c(7,8,9,10,11)]
        #test2pr <- preds[preds[,1] %in% c(unlist(folds[f,],use.names = FALSE)),c(7,8,9,10,11)]
        #test2n <- preds[preds[,1] %in% c(unlist(folds[f,],use.names = FALSE)),c(1)]
        test2pr <- preds[preds[,1] %in% c(folds[f,s]),c(7,8,9,10,11)]
        test2n <- preds[preds[,1] %in% c(folds[f,s]),c(1)]
        #colnames(train2pr) <- genes
        colnames(test2pr) <- genes
        #train2prc <- asdf[!preds[,1] %in% c(unlist(folds[f,],use.names = FALSE))]
        #test2prc <- asdf[preds[,1] %in% c(unlist(folds[f,],use.names = FALSE))]
        test2prc <- asdf[preds[,1] %in% c(folds[f,s])]
        #dtrain2pr = xgb.DMatrix(label=train2prc, data = as.matrix(train2pr))
        dtest2pr = xgb.DMatrix(label=test2prc, data = as.matrix(test2pr))
        #pr=round(predict(mdl,dtrain2pr))
        #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(train2prc, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
        #write.table(aa$table, paste('mil_cv/',wrun,'/trMIL.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
        pr=round(predict(mdl,dtest2pr))
        #teMIL[sss]=pr
        write.table(t(pr),paste('on_predictions/',prefix,'.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
        #aa=confusionMatrix(factor(pr,levels = c(0,1)), factor(test2prc, levels = c(0,1)), dnn = c("Pr","GT"), positive = "1")
        #write.table(aa$table, paste('on_predictions/',prefix,'_','teMIL.csv', sep = ''), sep = ",", append = T, col.names = FALSE, row.names = FALSE)
        
        # Write this
        #write.table(t(matrix(test2n,nrow=1)[1,]), paste('on_predictions/',prefix,'_','this.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
        #write.table(t(colSums(matrix(test2prc,nrow=1))), paste('on_predictions/',prefix,'_','this.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
        #write.table(t(colSums(matrix(pr,nrow=1))), paste('on_predictions/',prefix,'_','this.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
        
        
      }
    }
    sss=sss+1
  }
}
write.table(t(ns),paste('on_predictions/gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
write.table(t(gt),paste('on_predictions/gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
write.table(t(te),paste('on_predictions/gt.csv', sep = ''), sep = ',', append = T, col.names = FALSE, row.names = FALSE)
parallelStop()