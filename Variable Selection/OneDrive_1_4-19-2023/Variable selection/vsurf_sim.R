vsurf_simulation <- function(train.data, test.data, beta)
{  
  
  start.time = Sys.time()
  
  
  # Fit Sure Variable selecting using random forest
  X.train = data.matrix(X.train)
  vsurf_threshold <- VSURF_thres(X.train,factor(Y.train),verbose = T ,parallel = F)  
  vsurf_interpretation <- VSURF_interp(X.train,factor(Y.train),verbose = T,vars = vsurf_threshold$varselect.thres)
  imp_int = vsurf_interpretation$varselect.interp 
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
  miscla <- 1 - (sum(diag(cm))/sum(cm))
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
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  
  return(errors)
}
