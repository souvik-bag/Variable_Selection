library(stringr)
library(caret)


backward_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  
  # Start the time
  start.time <- Sys.time()
  # 
  # Fit the  model on the data
 glm_model <- glm(Y~., family = "binomial", data = train.data)
  
 # Fit forward variable selection
  
  backward_fit <- stepAIC(glm_model, direction = "backward")
  
  ## Selecting Important variables
  W <- as.matrix(coefficients(backward_fit))
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  ####Selected variables
  selected <- length(keep_X)
  ##select variables from weighted variables
  R = as.numeric(str_extract(keep_X, "\\d+"))
  selected1 = length(intersect(R,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(R,which(beta == 0)))/sum(beta == 0)
  
  # Estimated coefficients
  betahat <- rep(0,length(beta))
  betahat[R] <- coef(backward_fit)[-1]
  mse_beta <- mean((betahat - beta)^2)
  mae_beta <- mean(abs(betahat - beta))
  
  ## Prediction and performance of model with test data
  
  predicted_probs <- predict(backward_fit, newdata = as.data.frame(X.test), type = 'response')
  predicted_classes<- ifelse(predicted_probs > 0.5 , 1, 0)
  pred1 <- factor(predicted_classes)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  brier <- mean((Y.test - predicted_probs)^2)
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
                       time = as.numeric(time.taken), brier)
  return(errors)
}