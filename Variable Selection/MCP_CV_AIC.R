mcp_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  
  
  
  #minimum value of lambda 
  cv.mcp <- cv.ncvreg(X.train, Y.train, family = "binomial", penalty = "MCP")
  
  lambda_vec <- cv.mcp$lambda
  
  # Start time 
  start.time <- Sys.time()
  # Fit the  model on the training data
  mcp_model_cv <- ncvreg(X.train, Y.train, family = "binomial",
                      lambda = lambda_vec, penalty = "MCP")
  
  # Find the best lambda using AIC
  best_lambda=lambda_vec[which.min(AIC(mcp_model_cv))]
  # Now fit the model
  mcp_model <- ncvreg(X.train, Y.train, family = "binomial",
                      lambda = best_lambda, penalty = "MCP")
  ##########Important variables from MCP#######
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

