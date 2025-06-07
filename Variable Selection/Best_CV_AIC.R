best_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  start.time <- Sys.time()
  
  # Fit L0Learn by  to get the sequence of gamma and lambda 
  cvfit <- L0Learn.cvfit(X.train, Y.train, nFolds=5, loss = 'Logistic', penalty="L0L2", maxSuppSize=20, nGamma =10,
                         gammaMin=0.0001, gammaMax = 10)
  # get the list of lambda nad gamma
  lambda_list <- print(cvfit)
  
  # calculate the prediction probabilities
  
  for ( p in nrow(lambda_list)){
    predicted_probs[p] = predict(cvfit, newx = X.test, lambda = lambda_list[p,]$lambda, gamma =  lambda_list[p,]$gamma,  type = 'response')
    
  }
  
  
 
  cvfit3 <- L0Learn.fit(X.train, Y.train, loss = 'Logistic', penalty="L0L2", maxSuppSize=20, 
                        lambdaGrid =  list(0.175919)) 
  gamma_vec <- cvfit$gamma
  lambda_vec <- unlist(cvfit$lambda)
  AIC_values 
  
  # Compute the log-likelihood
  yhat <- predict(cvfit, newx = X.test)
  p <- 1 / (1 + exp(-yhat))
  loglik <- sum(Y.test * log(p[1,]) + (1 - Y.test) * log(1 - p[1,]))
  
  dbinom(as.numeric(yhat), prob=as.numeric(p), size=1, log=TRUE)
  
  
  L0Learn.fit(X.train, Y.train, loss = 'Logistic', penalty="L0L2", maxSuppSize=20, 
              lambdaGrid = lambda_vec)
  coef(cvfit,  lambda=2.45513e-02, gamma=0.0001)
  
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