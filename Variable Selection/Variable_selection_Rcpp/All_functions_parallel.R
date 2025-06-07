
### A-LASSO


alasso_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  # Make clusters
  
  registerDoParallel(detectCores()-1)
  
  #minimum value of lambda
  cv.ridge <- cv.glmnet(X.train, Y.train, alpha = 0, family = "binomial",  parallel = TRUE)
  # Fit the  model on the data
  ridge_model <- glmnet(X.train, Y.train, alpha = 0, family = "binomial",
                        lambda = cv.ridge$lambda.min, parallel = TRUE)
  betahat <- coef(ridge_model)[-1]
  best_ridge_weight= 1/ abs(as.matrix(betahat))
  # alasso_fit=glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight)
  alasso.cv=cv.glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight,type.measure = "deviance",nfold=10, parallel = TRUE)
  best_alasso=glmnet(X.train,Y.train,family="binomial",alpha=1,penalty.factor = best_ridge_weight,lambda=alasso.cv$lambda.min, intercept = F)
  ## Selecting Important variables
  W <- as.matrix(coef(best_alasso))
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ####Selected variables
  selected <- length(keep_X)
  ##select variables from weighted variables
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  
  
  betahat <- coef(best_alasso)[-1]
  mse_beta <- mean((betahat - beta)^2)
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  prediction <- best_alasso %>% predict(X.test, type = 'response')
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean((Y.test - prediction)^2)
  auc <- auc((Y.test), as.numeric(as.character(pred1)))
  ##error frame
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  
  errors <- data.frame(selected,
                       selected_Imp = round(selected1*100, 3),
                       selected_unImp = round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,3),
                       Recall = round( recall,3),
                       time = as.numeric(time.taken), brier, AUC = auc)
  return(errors)
}


best_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
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
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  
  betahat <- coeff[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  prediction <-  predict(cvfit, newx = X.test, lambda = optimalLambda, gamma =  optimalGamma,  type = 'response')
  
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean((Y.test - prediction[,1])^2)
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp = round(selected1*100,3),
                       selected_unImp = round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                      Accuracy = (accuracy*100),
                       Precision = round(precision,3),
                       Recall = round( recall,3),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
}


boruta_simulation <- function(X.train,Y.train,X.test,Y.test,beta){ 
  
  start.time <- Sys.time()
  
  # Fit Sure Variable selecting using random forest
  X.train = data.matrix(X.train)
  boruta <- Boruta(X.train,factor(Y.train), maxRuns = 1000)  
  
  a = parse_number(as.character(getSelectedAttributes(boruta)))
  selected <- length(getSelectedAttributes(boruta))
  selected1 = length(intersect(a,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(a,which(beta == 0)))/sum(beta == 0)
  # prediction using selected features
  varimp <-as.matrix( X.train[,a])
  rfimp <- randomForest(varimp,factor(Y.train))
  pred1 <- predict(rfimp,X.test[,a])
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, X.test[,a], na.action="na.impute", type = "prob")
  brier <- mean((Y.test - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp = round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = NA,
                       MAE = NA,
                       Accuracy = round(accuracy*100,3),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
}


elastic_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  # Create the data matrix
  train.data <- cbind(Y.train,X.train)
  train.data <- as.data.frame(train.data)
  train.data$Y.train <- as.factor(train.data$Y.train)
  # Create a cluster of worker processes
  cl <- makeCluster(10)
  
  # Register the cluster
  registerDoParallel(cl)
  
  # Define the cross-validation settings
  nfolds <- 10
  alpha_vals <- seq(0,1,length = 10)
  lambda_vals <-  seq(0.0001,1,length = 5)
  grid = expand.grid(alpha = alpha_vals,lambda = lambda_vals)
  
  # custom <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = F, allowParallel = F)
  # enet_model <- train(Y.train~.,train.data, method = 'glmnet', tuneGrid = expand.grid( alpha = seq(0,1,length = 10),lambda = seq(0.0001,1,length = 5)),
                      # trControl = custom)
  
  
 
  search <- foreach(i = alpha_vals, .combine = rbind) %dopar% {
    cv <- cv.glmnet(X.train, Y.train, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  # md3 <- glmnet(X.train, Y.train, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)

 elastic_model <- glmnet(X.train, Y.train, alpha =cv3$alpha, family = "binomial",
                          lambda = cv3$lambda.1se)
  
  # Estimated coefficients
  W <- as.matrix(coef(elastic_model))
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  selected <- length(keep_X)
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  
  betahat <- W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  ## Prediction and performance of model with test data
  
  prediction <- elastic_model %>% predict(X.test, type = 'response')
  
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  brier <- mean((Y.test - prediction)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  auc <- auc((Y.test), as.numeric(as.character(pred1)))
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp = round(selected1*100,3),
                       selected_unImp = round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                      Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier, AUC = auc)
  return(errors)
  
  
}


isis_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  #  Fit Sure Independence Screening
  
  sis_model <- SIS(X.train,Y.train,family = "binomial")
  
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
  
  R = W[-1]
  imp = row.names(W)[-1]
  
  betahat <- numeric(length(beta))
  betahat[unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))] = W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ##select variables from weighted variables
  selected1 = length(intersect(which(betahat != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(betahat != 0),which(beta == 0)))/sum(beta == 0)
  
  # prediction 
  prediction2 <- 1/ (1 + exp(-(intercept + (X.test%*%betahat))))
  pred1<- ifelse(prediction2 > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - prediction2)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}


lasso_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
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
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  ##select variables from weighted variables
  
  
  betahat <- W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test dpred1<- ifelse(prediction > 0.5 , 1, 0)
  prediction <- predict(lasso_model, X.test, type = "response")
  # prob = 1/ (1 + exp(-(coeff[1] + X.test%*% betahat)))
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - prediction)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp = round(selected1 * 100,3),
                       selected_unImp = round( selected2 * 100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                      Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}


mcp_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
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
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  ##select variables from weighted variables
  
  
  betahat <- W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  prediction <- mcp_model %>% predict(X.test, type = 'response')
  # probs <- 1/ (1 + exp(-(mcp_model$beta[1] + X.test %*% betahat)))
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - prediction)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  return(errors)
}


nta_simulation <- function(X.train,Y.train,X.test,Y.test,beta){ 
  
  start.time <- Sys.time()
  rf <- randomForest(X.train, factor(Y.train), importance = T, keep.inbag = T)
  
  vari = compVarImp(X.train, factor(Y.train),rf) 
  
  cvpvi <-   CVPVI(X.train, factor(Y.train), ntree = 500, k = 10, nPerm = 5, parallel = F)
  pval <- NTA(cvpvi$cv_varim)
  
  selected <- which (pval$pvalue <= 0.05)
  
  a <- length(selected)
  selected1 <- length(intersect(selected,which(beta != 0)))/sum(beta != 0)
  selected2 <- length(intersect(selected,which(beta == 0)))/sum(beta == 0)
  
  # prediction using selected features
  varimp <-as.matrix( X.train[,selected])
  rfimp <- randomForest(varimp,factor(Y.train))
  # cm <- table(Predicted = rfimp$predicted, Actual = Y.test)
  pred1 <- predict(rfimp,X.test[,selected])
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, X.test[,selected], na.action="na.impute", type = "prob")
  brier <- mean((Y.test - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected = a,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = NA,
                       MAE = NA,
                       Accuracy = round(accuracy*100,3),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
  
}


pimp_simulation <- function(X.train,Y.train,X.test,Y.test,beta){ 
  
  start.time <- Sys.time()

  
  # Register parallel backend
  cl <- makeCluster( detectCores() -1 )
  cores = detectCores() -1 
  registerDoParallel(cl)
  
  # Train Random Forest model in parallel
  rf_model <- system.time( foreach(ntree = rep(500, cores), .combine = combine, .packages = "randomForest") %dopar% {
    randomForest(X.train, factor(Y.train), ntree = ntree, importance = T, keep.inbag = T, mtry = 3)
  })
  
  # Stop parallel backend
  stopCluster(cl)
  
  
  
  
  
  
  rf <- system.time(randomForest(X.train, factor(Y.train), importance = T, keep.inbag = T))
  
  pimp <- PIMP(X.train, factor(Y.train),rf)
  
  pimptest <- PimpTest(pimp, para = F)
  
  selected <- which(pimptest$pvalue <= 0.05)
  a <- length(selected)
  selected1 <- length(intersect(selected,which(beta != 0)))/sum(beta != 0)
  selected2 <- length(intersect(selected,which(beta == 0)))/sum(beta == 0)
  
  # prediction using selected features
  varimp <-as.matrix( X.train[,selected])
  rfimp <- randomForest(varimp,factor(Y.train))
  pred1 <- predict(rfimp, X.test[,selected])
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy =  (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  prediction = predict(rfimp, X.test[,selected], na.action="na.impute", type = "prob")
  brier <- mean((Y.test - prediction[,2])^2)
  ##error frame
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(Selected = a,
                       Selected_Imp = round(selected1*100,3),
                       Selected_unImp =  round(selected2*100,3),
                       MSE = NA,
                       MAE = NA,
                       Accuracy = round(accuracy*100 , 3),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  return(errors)
  
}


rrf_simulation <- function(X.train,Y.train,X.test,Y.test,beta){ 
  
  
  start.time <- Sys.time()
  
  #ordinary random forest.
  rf <- RRF(X.train,factor(Y.train), flagReg = 0, importance = TRUE)
  impRF <- rf$importance
  impRF <- impRF[,"MeanDecreaseGini"]
  #guided regularized random forest
  imp <- impRF/(max(impRF))#normalize the importance score
  gamma <- 0.5
  coefReg <- (1-gamma)+gamma*imp #weighted average
  grrf <- RRF(X.train,factor(Y.train),coefReg=coefReg,  flagReg=1, importance = TRUE)
  
  # Selected set of features used for checking prediction accuracy
  a <- grrf$feaSet
  selected <- length(a)
  selected1 = length(intersect(a,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(a,which(beta == 0)))/sum(beta == 0)
  # when Not regularized this may happen
  
  # prediction using selected features
  varimp <- as.data.frame(X.train[,a])
  rfimp <- randomForest(x = varimp, y = factor(Y.train))
  
  pred1 <- predict(rfimp,X.test[,a])
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  prediction = predict(rfimp, X.test[,a], na.action="na.impute", type = "prob")
  brier <- mean((Y.test - prediction[,2])^2)
  
  ##error frame
  
  
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = NA,
                       MAE = NA,
                       Accuracy = round(accuracy*100,3),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}


scad_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
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
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  ##select variables from weighted variables
  # selected1 <- sum(ifelse(as.numeric(substr(keep_X,2,5)) <= ,1,0))
  
  betahat <- W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  prediction <- scad_model %>% predict(X.test, type = 'response')
  
  pred1 <- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - prediction)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  ##error frame
  
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  
  return(errors)
  
}


sis_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  #  Fit Sure Independence Screening
  
  sis_model <- SIS(X.train,Y.train,family = "binomial", varISIS = "vanilla", penalty = "SCAD", iter = F)
  
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
  
  R = W[-1]
  imp = row.names(W)[-1]
  
  betahat <- numeric(length(beta))
  betahat[unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))] = W[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ##select variables from weighted variables
  selected1 = length(intersect(which(betahat != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(betahat != 0),which(beta == 0)))/sum(beta == 0)
  
  # prediction 
  prediction2 <- 1/ (1 + exp(-(intercept + (X.test%*%betahat))))
  pred1<- ifelse(prediction2 > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - prediction2)^2)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
  
}


sivs_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  sivs_model <- sivs(X.train,factor(Y.train), family = "binomial",  return.fits = T, verbose = "general", progressbar = T,parallel.cores = "grace" ) 
  imp = sivs::suggest(sivs_model, strictness = 0.5)
  
  ## glmnet using sivs
  
  sivs_glmnet_model <- glmnet::cv.glmnet(x = data.matrix(X.train[,imp]),
                                         y = factor(Y.train),
                                         family = "binomial")
  coeff <- coef(sivs_glmnet_model)
  
  W <- as.matrix(coeff)
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ###Selected variables
  selected <- length(keep_X)
  R = W[-1]
  imp_position = unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))
  selected1 = length(intersect(imp_position,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(imp_position,which(beta == 0)))/sum(beta == 0)  
  
  betahat <- numeric(length(beta))
  betahat[unlist( sapply(imp, FUN = function(xx) which(colnames(X.train) == xx)))] = R
  
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  prediction <- predict(object = sivs_glmnet_model,
                        newx = X.test[,imp],
                        s = "lambda.min",
                        type = "response")
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean((Y.test - prediction)^2)
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp = round(selected1*100,3),
                       selected_unImp = round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
  
}


sparse_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  start.time <- Sys.time()
  
  
  # Fit L0Learn by Cross Validation
  
  cvfit <- L0Learn.cvfit(X.train, Y.train, nFolds=5, loss = 'Logistic', penalty="L0", maxSuppSize=20, nGamma =10,
                         gammaMin=0.0001, gammaMax = 10)
  
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
  ##select variables from weighted variables
  R = W[-1]
  selected1 = length(intersect(which(R != 0),which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(which(R != 0),which(beta == 0)))/sum(beta == 0)
  
  betahat <- coeff[-1]
  mse_beta <- (mean((betahat - beta)^2))
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  prediction <-  predict(cvfit, newx = X.test, lambda = optimalLambda, gamma =  optimalGamma,  type = 'response')
  
  pred1<- ifelse(prediction > 0.5 , 1, 0)
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean((Y.test - prediction[,1])^2)
  ##error frame
  # Time taken
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp = round(selected2*100,3),
                       MSE = round(mse_beta,2),
                       MAE = round(mae_beta,2),
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken),brier)
  
  return(errors) 
}


vsurf_simulation <- function(X.train,Y.train,X.test,Y.test,beta){  
  
  start.time = Sys.time()
  
  
  # Fit Sure Variable selecting using random forest
  X.train = data.matrix(X.train)
  vsurf_model <- VSURF(X.train,factor(Y.train),verbose = F ,parallel = T, ncores = detectCores()-1)  
  imp_int = vsurf_model$varselect.interp 
  forest_imp <- randomForest(X.train[,imp_int], factor(Y.train), ntree = 500)
  
  prediction =  predict(forest_imp,X.test, type = "prob")
  
  pred <- prediction[,2]
  # pred1<- ifelse(pred > 0.5 , 1, 0)
  pred1 <-  predict(forest_imp,X.test, type = "response")
  
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - pred)^2)
  #  cm <- table(Predicted = pred_t, Actual = Y.test)
  # # 
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  #selected
  
  selected = length(imp_int)
  selected1 = length(intersect(imp_int,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(imp_int,which(beta == 0)))/sum(beta == 0)
  # Time taken
  end.time <- Sys.time() 
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       selected_Imp =  round(selected1*100,3),
                       selected_unImp =  round(selected2*100,3),
                       MSE = NA,
                       MAE = NA,
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  
  return(errors)
}
