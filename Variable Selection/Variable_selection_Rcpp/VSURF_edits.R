install.packages("VSURF")
library(VSURF)

sourceCpp("VSURFcpp.cpp")


# set.seed(rpts)
X = matrix(rnorm(setups[i,]$n*setups[i,]$m,0,1),nrow = setups[i,]$n, ncol = setups[i,]$m)

beta <- numeric(setups[i,]$m)
beta[sample(setups[i,]$m,setups[i,]$k)] = rnorm(setups[i,]$k,0,3)
beta <- ifelse((beta>0 & beta <0.5),beta + 0.5, ifelse((beta <0 & beta > -0.5), beta-0.5,beta))

####we get the value of log(p/(1-p)) and errors are taken normal
logitp = X%*%beta + rnorm(setups[i,]$n)
######by some manipulation we get the value of p(y=1)
p = exp(logitp)/(1 +exp(logitp) )
#Using this value of p, we get the binary response variable
Y <- rbinom(setups[i,]$n,1,p)

#####merge X and Y to create one dataframe
sim_data <- as.data.frame(cbind(Y,X))
# set_names(sim_data[,-1],paste0("X",1:10))

names = c("Y",paste0("X",1:setups[i,]$m))
colnames(sim_data) <- names
## Train Test split
training.sample <- sim_data$Y %>% createDataPartition(p = 0.8 , list = FALSE) 

train.data <- sim_data[training.sample,]
test.data <- sim_data[-training.sample,]

#######Now we fit lasso
Y.train <- train.data$Y
Y.test <- test.data$Y
#model matrix
X.train <- (model.matrix(Y.train~.,data = train.data))[,-c(1,2)]
X.test <- model.matrix(Y.test~.,data = test.data)[,-c(1,2)]


## VSURF

X.train = data.matrix(X.train)
vsurf_model <- VSURF(X.train,factor(Y.train),verbose = F ,parallel = F) 
vsurf_model_cpp <- VSURFcpp(x = X.train, y = factor(Y.train), ntree = 500, verbose = FALSE, parallel= FALSE)
imp_int = vsurf_model$varselect.pred
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
accuracy <- (cm[1,1]+cm[2,2])/ (cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1])
precision <- cm[1,1]/(cm[1,1]+cm[1,2])
recall <- cm[1,1]/(cm[1,1]+cm[2,1])
auc <- auc((Y.test), as.numeric(as.character(pred1)))
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
                     Accuracy = (accuracy*100),
                     Precision = round(precision,2),
                     Recall = round( recall,2),
                     time = as.numeric(time.taken), brier,
                     AUC = auc)

