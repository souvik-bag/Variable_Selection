library(sda)
library(dplyr)
library(glmnet)
library(tidyverse)
library(sivs)
library(L0Learn)
library(ncvreg)
library(caret)
library(SIS)
library(Boruta)
library(VSURF)
library(vita)
library(RRF)
library(randomForest)
data(singh2002)

dim(singh2002$x)
length(singh2002$y)
singh2002$y <- recode_factor(singh2002$y, healthy = 0, 
                             cancer = 1)
singh2002$y
 X = singh2002$x
 Y = singh2002$y


vsurf_simulation_r  <- function(X,Y){  
  
  start.time = Sys.time()
  
  
  # Fit Sure Variable selecting using random forest
  X = data.matrix(X)
  vsurf <- VSURF(X,factor(Y),verbose = T ,parallel = T, ncores = detectCores() - 1)
  imp_int = vsurf$varselect.interp 
  forest_imp <- randomForest(X[,imp_int], factor(Y), ntree = 500)
  
  prediction =  predict(forest_imp,X[,imp_int], type = "prob")
  
  pred <- prediction[,2]
  pred1 <-  predict(forest_imp,X[,imp_int], type = "response")
  
  cm <- confusionMatrix(data = pred1, reference = factor(Y))
  cm <- as.matrix(cm)
  brier <- mean((as.numeric(levels(Y))[Y] - pred)^2)
  miscla <- 1 - (sum(diag(cm))/sum(cm))
  precision <- cm[1,1]/(cm[1,1]+cm[1,2])
  recall <- cm[1,1]/(cm[1,1]+cm[2,1])
  
  #selected
  
  selected = length(imp_int)
  
  end.time <- Sys.time() 
  time.taken <- difftime(end.time,start.time,units = "min")
  ##error frame
  errors <- data.frame(selected,
                       Missclassification_rate = miscla,
                       Precision = round(precision,2),
                       Recall = round( recall,2),
                       time = as.numeric(time.taken), brier)
  
  
  
  
  return(errors)
}
