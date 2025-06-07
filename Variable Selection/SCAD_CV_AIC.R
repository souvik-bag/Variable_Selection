scad_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  #minimum value of lambda 
  cv.scad <- cv.ncvreg(X.train, Y.train, family = "binomial", penalty = "SCAD")
  # get the lambda values
  lambda_vec <- cv.mcp$lambda
  # Start time 
  start.time <- Sys.time()
  # Fit the  model on the training data
  scad_model_cv <- ncvreg(X.train, Y.train, family = "binomial",
                       lambda = lambda_vec, penalty = "SCAD")
  # Find the best lambda using AIC
  best_lambda=lambda_vec[which.min(AIC(scad_model_cv))]
  # Now fit the model
  scad_model <- ncvreg(X.train, Y.train, family = "binomial",
                      lambda = best_lambda, penalty = "SCAD")
  
  
  ##########Important variables from SCAD#######
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
