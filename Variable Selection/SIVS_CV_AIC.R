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
