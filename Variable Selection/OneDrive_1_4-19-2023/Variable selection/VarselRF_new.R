varselrf_simulation <- function(X.train,Y.train,X.test,Y.test,beta){
  
  
  # Start the time
  start.time <- Sys.time()
  
  ## Fit a VarSelRF model
  model <-varSelRF(xdata = X.train, Class = as.factor(Y.train), ntree = 500, ntreeIterat = 500)
  
  ###Selected variables
  selected_vars <- model$selected.vars
  selected <-  length(as.numeric(sub("X", "", selected_vars)))
  
  # Selection metrics
  selected1 = length(intersect(selected,which(beta != 0)))/sum(beta != 0)
  selected2 = length(intersect(selected,which(beta == 0)))/sum(beta == 0)
  
  ## Prediction and performance of model 
fitted_model <- randomForest(y=as.factor(Y.train), x = subset(X.train,select=model$selected.vars), ntree=model$ntreeIterat, xtest=subset(X.test,select=model$selected.vars))
  
  pred1<- fitted_model$predicted
  pred1 <- factor(pred1)
  cm <- confusionMatrix(data = pred1, reference = factor(Y.test))
  cm <- as.matrix(cm)
  brier <- mean((Y.test - fitted_model$votes[,2])^2)
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
                       Accuracy = (accuracy*100),
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  return(errors)
}


varselrf_simulation(X.train,Y.train,X.test,Y.test,beta)
