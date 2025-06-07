
library(sda)
library(dplyr)
library(glmnet)
library(tidyverse)
library(sivs)
library(L0Learn)
library(ncvreg)
library(caret)
library(SIS)
library(Boruta)
library(VSURF)
library(vita)
library(RRF)
library(randomForest)
data(singh2002)

dim(singh2002$x)
length(singh2002$y)
singh2002$y <- recode_factor(singh2002$y, healthy = 0, 
                                          cancer = 1)

rpts = 1
rptns = 20

real_lasso = matrix(nrow = rptns,ncol = 6)
real_enet = matrix(nrow = rptns,ncol = 6)
real_alasso = matrix(nrow = rptns, ncol = 6)
real_sivs = matrix(nrow = rptns,ncol = 6)
real_sparse = matrix(nrow = rptns,ncol = 6)
real_best = matrix(nrow = rptns,ncol = 6)
real_scad = matrix(nrow = rptns,ncol = 6)
real_mcp = matrix(nrow = rptns,ncol = 6)
real_sis = matrix(nrow = rptns,ncol = 6)
real_isis = matrix(nrow = rptns,ncol = 6)
real_boruta = matrix(nrow = rptns,ncol = 6)
real_vsurf = matrix(nrow = rptns,ncol = 6)
real_rrf = matrix(nrow = rptns,ncol = 6)
real_pimp =  matrix(nrow = rptns,ncol = 6)
real_nta =  matrix(nrow = rptns,ncol = 6)


for(rpts in 1:rptns){

## Train Test split
training.sample <-  singh2002$y %>% createDataPartition(p = 0.8 , list = FALSE) 

X.train = singh2002$x[training.sample,] 
Y.train = singh2002$y[training.sample] 
X.test = singh2002$x[-training.sample,] 
Y.test = singh2002$y[-training.sample] 


## 1 - lasso ##
lasso_simulation_r <-function(X.train,Y.train, X.test,Y.test){
  
  start.time <- Sys.time()
  
  #minimum value of lambda 
  cv.lasso <- cv.glmnet(X.train, Y.train, alpha = 1, family = "binomial", intercept = F)
  # Fit the  model on the training data
  lasso_model <- glmnet(X.train, Y.train, alpha = 1, family = "binomial",
                        lambda = cv.lasso$lambda.min, intercept = F)
  coeff <- coef(lasso_model)
  ## Selecting Important variables
  
  W <- as.matrix(coeff)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ###Selected variables
  selected <- length(keep_X)
  
  ## Prediction and performance of model with test dpred1<- ifelse(prediction > 0.5 , 1, 0)
  prediction <- predict(lasso_model, X.test, type = "response")
  # prob = 1/ (1 + exp(-(coeff[1] + X.test%*% betahat)))
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}

alasso_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  start.time <- Sys.time()
  
  #minimum value of lambda
  cv.ridge <- cv.glmnet(X.train, Y.train, alpha = 0, family = "binomial")
  # Fit the  model on the data
  ridge_model <- glmnet(X.train, Y.train, alpha = 0, family = "binomial",
                        lambda = cv.ridge$lambda.min)
  betahat <- coef(ridge_model)[-1]
  best_ridge_weight= 1/ abs(as.matrix(betahat))
  alasso_fit=glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight)
  alasso.cv=cv.glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight,type.measure = "deviance",nfold=10)
  best_alasso=glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight,lambda=alasso.cv$lambda.min, intercept = F)
  ## Selecting Important variables
  W <- as.matrix(coef(best_alasso))
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ####Selected variables
  selected <- length(keep_X)
  
  ## Prediction and performance of model with test data
  
  prediction <- best_alasso %>% predict(X.test, type = 'response')
  
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  ##error frame
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  return(errors)
}


best_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  
  start.time <- Sys.time()
  
  # Fit L0Learn by Cross Validation
  
  cvfit <- L0Learn.cvfit(X.train, Y.train, nFolds=5, loss = 'Logistic', penalty="L0L2", maxSuppSize=20, nGamma =10,
                         gammaMin=0.0001, gammaMax = 10)
  
  cvErrors <- lapply(cvfit$cvMeans, min)
  
  index <- which.min(cvErrors)
  
  optimalGammaIndex = index # index of the optimal gamma identified previously
  optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])
  optimalLambda = cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
  optimalGamma = cvfit$fit$gamma[index] 
  
  # Coefficients of fitted model
  
  coeff <- coef(cvfit, lambda=optimalLambda, gamma= optimalGamma)
  
  ##########Important variables from Lasso#######
  W <- as.matrix(coeff)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "Intercept"]
  
  ####Selected variables
  selected <- length(keep_X)
  ##select variables from weighted variables
 
  
  ## Prediction and performance of model with test data
  
  prediction <-  predict(cvfit, newx = X.test, lambda = optimalLambda, gamma =  optimalGamma,  type = 'response')
  
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean(((as.numeric(levels(Y.test))[Y.test])- prediction)^2)
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
}


elastic_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  start.time <- Sys.time()
   
  train.data$Y<- as.factor(train.data$Y)
  
  
  train.data$Y<- as.factor(train.data$Y)
  custom <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = F, allowParallel = F)
  enet_model <- train(Y~.,train.data, method = 'glmnet', tuneGrid = expand.grid( alpha = seq(0,1,length = 10),lambda = seq(0.0001,1,length = 10)),
                      trControl = custom)
  elastic_model <- glmnet(X.train, Y.train, alpha = enet_model$bestTune$alpha, family = "binomial",
                          lambda = enet_model$bestTune$lambda)
  
  
  W <- as.matrix(coef(elastic_model))
  
  
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  selected <- length(keep_X)
 
  ## Prediction and performance of model with test data
  
  prediction <- elastic_model %>% predict(X.test, type = 'response')
  
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  return(errors)
  
  
}

sparse_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  
  start.time <- Sys.time()
  
  
  # Fit L0Learn by Cross Validation
  
  cvfit <- L0Learn.cvfit(X.train, Y.train, nFolds=10, loss = 'Logistic', penalty="L0", maxSuppSize=20, nGamma =100,
                         gammaMin=0.001, gammaMax = 10)
  
  cvErrors<- lapply(cvfit$cvMeans, min)
  
  index <- which.min(cvErrors)
  
  optimalGammaIndex = index # index of the optimal gamma identified previously
  optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])
  optimalLambda = cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
  optimalGamma = cvfit$fit$gamma[index] 
  
  # Coefficients of fitted model
  
  coeff <- coef(cvfit, lambda=optimalLambda, gamma= optimalGamma)
  
  ##########Important variables from Lasso#######
  W <- as.matrix(coeff)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "Intercept"]
  
  ####Selected variables
  selected <- length(keep_X)

  ## Prediction and performance of model with test data
  
  prediction <-  predict(cvfit, newx = X.test, lambda = optimalLambda, gamma =  optimalGamma,  type = 'response')
  
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  ##error frame
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
}

sivs_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  start.time <- Sys.time()
  
  sivs_model <- sivs(X.train,factor(Y.train), family = "binomial",  return.fits = T, verbose = "general", progressbar =  T
                     , iter.count = 50, nfolds = 5 ,parallel.cores = NULL) 
  imp = sivs::suggest(sivs_model, strictness = 0.5)
  imp <- as.character(imp)
  imp <- parse_number(imp)
  ## glmnet using sivs
  
  sivs_glmnet_model <- glmnet::cv.glmnet(x = X.train[,imp],
                                         y = factor(Y.train),
                                         family = "binomial")
  coeff <- coef(sivs_glmnet_model)
  
  W <- as.matrix(coeff)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ###Selected variables
  selected <- length(keep_X)
 
  
  prediction <- predict(object = sivs_glmnet_model,
                        newx = X.test[,imp],
                        s = "lambda.min",
                        type = "response")
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
  
}

mcp_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  
  start.time <- Sys.time()
  
  
  #minimum value of lambda 
  cv.mcp <- cv.ncvreg(X.train, Y.train, family = "binomial", penalty = "MCP")
  best_lambda=cv.mcp$lambda.min
  
  # Fit the  model on the training data
  mcp_model <- ncvreg(X.train, Y.train, family = "binomial",
                      lambda = best_lambda, penalty = "MCP")
  
  ##########Important variables from Lasso#######
  W <- as.matrix(mcp_model$beta)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  
  ####Selected variables
  selected <- length(keep_X)
  
  ## Prediction and performance of model with test data
  
  prediction <- mcp_model %>% predict(X.test, type = 'response')
  # probs <- 1/ (1 + exp(-(mcp_model$beta[1] + X.test %*% betahat)))
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  return(errors)
}

scad_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  
  start.time <- Sys.time()
  
  
  #minimum value of lambda 
  cv.scad <- cv.ncvreg(X.train, Y.train, family = "binomial", penalty = "SCAD")
  
  best_lambda=cv.scad$lambda.min
  
  # Fit the  model on the training data
  scad_model <- ncvreg(X.train, Y.train, family = "binomial",
                       lambda = best_lambda, penalty = "SCAD")
  
  ##########Important variables from Lasso#######
  W <- as.matrix(scad_model$beta)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  
  ####Selected variables
  selected <- length(keep_X)
  
  ## Prediction and performance of model with test data
  
  prediction <- scad_model %>% predict(X.test, type = 'response')
  
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  
  return(errors)
  
}

sis_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  start.time <- Sys.time()
  Y.train_num = as.numeric(as.character(Y.train))
  #  Fit Sure Independence Screening
  
  sis_model <- SIS(X.train,Y.train_num,family = "binomial", iter = F)
  
  ## Selecting Important variables
  #sis_model
  intercept <- as.matrix(sis_model$coef.est)[1,]
  intercept = as.numeric(intercept)
  W <-as.matrix( sis_model$coef.est)
  keep_X <- rownames(W)
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ####Selected variables
  selected <- length(keep_X)
  ##select variables from weighted variables
  imp = row.names(W)[-1]
  imp_position = parse_number(as.character(imp))
  betahat <- numeric(ncol(X.train))
  # betahat[unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))] = W[-1]
 betahat[imp_position]=W[-1]
  
  # prediction 
  prediction2 <- 1/ (1 + exp(-(intercept + (X.test%*%betahat))))
  pred1<- ifelse(prediction2 > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction2)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
  
}

isis_simulation_r <- function(X.train,Y.train, X.test,Y.test){
  start.time <- Sys.time()
  Y.train_num = as.numeric(as.character(Y.train))
  #  Fit Sure Independence Screening
  
  sis_model <- SIS(X.train,Y.train_num,family = "binomial", varISIS = "cons", penalty = "MCP", iter = T, iter.max = 10000)
  
  ## Selecting Important variables
  #sis_model
  intercept <- as.matrix(sis_model$coef.est)[1,]
  intercept = as.numeric(intercept)
  W <-as.matrix( sis_model$coef.est)
  keep_X <- rownames(W)
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ####Selected variables
  selected <- length(keep_X)
  ##select variables from weighted variables
  imp = row.names(W)[-1]
  imp_position = parse_number(as.character(imp))
  betahat <- numeric(ncol(X.train))
  # betahat[unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))] = W[-1]
  betahat[imp_position]=W[-1]
  
  
  # prediction 
  prediction2 <- 1/ (1 + exp(-(intercept + (X.test%*%betahat))))
  pred1<- ifelse(prediction2 > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean(((as.numeric(levels(Y.test))[Y.test]) - prediction2)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
  
  
  
}

boruta_simulation_r <- function(X.train,Y.train, X.test,Y.test){ 
  
  start.time <- Sys.time()
  
  
  # Fit Sure Variable selecting using random forest
  X.train = data.matrix(X.train)
  boruta <- Boruta(X.train,factor(Y.train), maxRuns = 1000)  
  # a <- as.numeric(substr(getSelectedAttributes(boruta),2,5))
  a = parse_number(as.character(getSelectedAttributes(boruta)))
  selected <- length(getSelectedAttributes(boruta))
  
  varimp <-as.matrix( X.test[,a])
  rfimp <- randomForest(varimp,factor(Y.test))
  pred1 <- rfimp$predicted
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))     
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, varimp, na.action="na.impute", type = "prob")
  brier <- mean((as.numeric(levels(Y.test))[Y.test] - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
}

vsurf_simulation_r  <- function(X.train,Y.train, X.test,Y.test){  
  
  start.time = Sys.time()
  
  
  # Fit Sure Variable selecting using random forest
  X.train = data.matrix(X.train)
  vsurf_threshold <- VSURF_thres(X.train,factor(Y.train),verbose = T ,parallel = T)  
  vsurf_interpretation <- VSURF_interp(X.train,factor(Y.train),verbose = T,vars = vsurf_threshold$varselect.thres,parallel = T )
  imp_int = vsurf_interpretation$varselect.interp 
  forest_imp <- randomForest(X.train[,imp_int], factor(Y.train), ntree = 500)
  
  prediction =  predict(forest_imp,X.test, type = "prob")
  
  pred <- prediction[,2]
  # pred1<- ifelse(pred > 0.5 , 1, 0)
  pred1 <-  predict(forest_imp,X.test, type = "response")
  
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((as.numeric(levels(Y.test))[Y.test] - pred)^2)
  #  cm <- table(Predicted = pred_t, Actual = Y.test)
  # # 
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  #selected
  
  selected = length(imp_int)
  
  end.time <- Sys.time() 
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  
  return(errors)
}

rrf_simulation_r <- function(X.train,Y.train, X.test,Y.test){ 
  
  
  start.time <- Sys.time()
  
  #ordinary random forest.
  rf <- RRF(X.train,factor(Y.train), flagReg = 0, importance = TRUE)
  impRF <- rf$importance
  impRF <- impRF[,"MeanDecreaseGini"]
  #guided regularized random forest
  imp <- impRF/(max(impRF))#normalize the importance score
  gamma <- 0.5
  coefReg <- (1-gamma)+gamma*imp #weighted average
  grrf <- RRF(X.test,factor(Y.test),coefReg=coefReg,  flagReg=1, importance = TRUE)
  
  # Selected set of features used for checking prediction accuracy
  a <- grrf$feaSet
  selected <- length(a)
  # when Not regularized this may happen
  # 
  # rf2 <-RRF(X.train,factor(Y.train), flagReg = 1, importance = TRUE)
  # rf2$feaSet
  
  # prediction using selected features
  varimp <- as.data.frame(X.test[,a])
  rfimp <- randomForest(x = varimp, y = factor(Y.test))
  
  pred1 <- rfimp$predicted
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  prediction = predict(rfimp, varimp, na.action="na.impute", type = "prob")
  brier <- mean((as.numeric(levels(Y.test))[Y.test]- prediction[,2])^2)
  
  ##error frame
  
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}

pimp_simulation_r <- function(X.train,Y.train, X.test,Y.test){ 
  
  start.time <- Sys.time()
  rf <- randomForest(X.train, factor(Y.train), importance = T, keep.inbag = T)
  
  pimp <- PIMP(X.train, factor(Y.train),rf)
  
  pimptest <- PimpTest(pimp, para = F)
  
  selected <- which(pimptest$pvalue <= 0.05)
  a <- length(selected)
  
  # prediction using selected features
  varimp <-as.matrix( X.test[,selected])
  rfimp <- randomForest(varimp,factor(Y.test))
  pred1 <- rfimp$predicted
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  miscla <- 1 - (sum(diag(cm))/sum(cm))     
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, varimp, na.action="na.impute", type = "prob")
  brier <- mean((as.numeric(levels(Y.test))[Y.test] - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected = a,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
  
}

nta_simulation_r <- function(X.train,Y.train, X.test,Y.test){ 
  
  start.time <- Sys.time()
  rf <- randomForest(X.train, factor(Y.train), importance = T, keep.inbag = T)
  
  vari = compVarImp(X.train, factor(Y.train),rf) 
  # pimp <- PIMP(X.train, factor(Y.train),rf)
  # 
  # pimptest <- PimpTest(pimp, para = F)
  # 
  # selected <- which(pimptest$pvalue <= 0.05)
  cvpvi <-   CVPVI(X.train, factor(Y.train), ntree = 500, k = 10, nPerm = 5)
  pval <- NTA(cvpvi$cv_varim)
  
  selected <- which (pval$pvalue <= 0.05)
  
  a <- length(selected)
 
  # prediction using selected features
  varimp <-as.matrix( X.test[,selected])
  rfimp <- randomForest(varimp,factor(Y.test))
  # cm <- table(Predicted = rfimp$predicted, Actual = Y.test)
  pred1 <- rfimp$predicted
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  miscla <- 1 - (sum(diag(cm))/sum(cm))     
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, varimp, na.action="na.impute", type = "prob")
  brier <- mean((as.numeric(levels(Y.test))[Y.test] - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected = a,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
  
}

real_lasso <- as.numeric(lasso_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_alasso <- as.numeric(alasso_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_best <- as.numeric(best_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_elestic <- as.numeric(elastic_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_sparse <- as.numeric(sparse_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_sivs <- as.numeric(sivs_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_mcp <- as.numeric(mcp_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_scad <- as.numeric(scad_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_sis <- as.numeric(sis_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_isis <- as.numeric(isis_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_boruta <- as.numeric(boruta_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_vsurf <- as.numeric(vsurf_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_rrf <- as.numeric(rrf_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_pimp <- as.numeric(pimp_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))
real_nta <- as.numeric(nta_simulation_r(X.train = X.train, Y.train = Y.train, X.test = X.test, Y.test = Y.test))

print("Number of Iteration")
print(rpts)
}
real_result <- rbind(real_alasso,real_elestic)
