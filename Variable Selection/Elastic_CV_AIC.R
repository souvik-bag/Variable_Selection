elastic_net_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  
  # Calculate AIC for glmnet object
  
  aic_glmnet <- function(fit){
    
    tLL <- - deviance(fit)
    k <- fit$df + 1
    n <- nobs(fit)
    AICc <- -tLL+2*k
    return(AICc)
  }
  
  # Start the time
  start.time <- Sys.time()
  # First we get the sequence of lambdas
  cv.lasso <- cv.glmnet(X.train, Y.train, alpha = 1, family = "binomial")
  
  # Run a for loop to calculate AIC for different lambdas
  
  lambda_values = seq(0.01,10, length.out = 1000) # If want to provide lambda sequence by ourselves
  alpha_values <- seq(0.01, 0.99, by = 0.5)
  
  # Create a grid of alpha and lambda values
  grid <- expand.grid(alpha = alpha_values, lambda = lambda_values)
  
  AIC_values = data.frame(alpha = NULL, lambda = NULL , AIC = NULL)
  
  for (l in lambda_values){
    for (a in alpha_values){
    aic.lasso <- glmnet(X.train, Y.train, alpha = a, family = "binomial", lambda = l)
    
    new_row = data.frame(alpha = a, lambda = l, AIC =   aic_glmnet(aic.lasso))
    AIC_values <- rbind(AIC_values, new_row)
    
  }}
  
  
  #minimum value of lambda
  lambda.best = AIC_values$lambda[which.min(AIC_values$AIC)]
  #minimum value of alpha
  alpha.best = AIC_values$alpha[which.min(AIC_values$AIC)]
 
   # Fit the  model on the data
  elastic_net_model <- glmnet(X.train, Y.train, alpha = alpha.best, family = "binomial",
                        lambda = lambda.best)
  coeff <- coef(elastic_net_model)
  
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

