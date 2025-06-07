sis_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  start.time <- Sys.time()
  
  #  Fit Sure Independence Screening
  
  sis_model <- SIS(X.train,Y.train,family = "binomial", varISIS = "vanilla", penalty = "SCAD", iter = F, tune = 'aic')
  
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
