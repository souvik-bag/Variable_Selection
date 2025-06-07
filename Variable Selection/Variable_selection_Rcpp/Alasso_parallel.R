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
