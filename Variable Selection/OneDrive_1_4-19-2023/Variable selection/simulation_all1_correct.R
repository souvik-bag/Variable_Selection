
## Uncorrelated Simulation, put i=4.


# LOAD LIBRARIES

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
library(parallel)
library(magrittr)
library(caret)

## create setup frame
n = c(100,100,100,10000)
m = c(10,100,1000,1000)
k = c(3,5,10,200)
setups = data.frame(n,m,k)
rptns = 2
rpts = 1


# Select the value of i

i <- 1

r_lasso14 = matrix(nrow = rptns,ncol = 10)
r_enet14 = matrix(nrow = rptns,ncol = 10)
r_alasso14 = matrix(nrow = rptns, ncol=10)
r_sivs14 = matrix(nrow = rptns,ncol = 10)
r_sparse14 = matrix(nrow = rptns,ncol = 10)
r_best14 = matrix(nrow = rptns,ncol = 10)
r_scad14 = matrix(nrow = rptns,ncol = 10)
r_mcp14 = matrix(nrow = rptns,ncol = 10)
r_sis14 = matrix(nrow = rptns,ncol = 10)
r_isis14 = matrix(nrow = rptns,ncol = 10)
r_boruta14 = matrix(nrow = rptns,ncol = 10)
r_vsurf14 = matrix(nrow = rptns,ncol = 10)
r_rrf14 = matrix(nrow = rptns,ncol = 10)
r_pimp14 =  matrix(nrow = rptns,ncol = 10)
r_nta14 =  matrix(nrow = rptns,ncol = 10)



for (rpts in 1:rptns)
{
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


  ## Start counting one iteration
  
  begin <- Sys.time()
  
  # ## Boruta ## Lasso
  expt = try(expr = {
    r_lasso14[rpts,] = as.numeric(lasso_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta ))
  },silent = T)
  ## Elastic net
  expt = try(expr = {
    r_enet14 [rpts,] = as.numeric(elastic_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  ## Awdaptive lasso
  expt = try(expr = {
    r_alasso14[rpts,] = as.numeric( alasso_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  
  # ## L0L2
  expt = try(expr = {
    r_best14[rpts,] = as.numeric(best_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## Sparse Step
  expt = try(expr = {
    r_sparse14[rpts,] = as.numeric(sparse_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## SCAD
  expt = try(expr = {
    r_scad14[rpts,] = as.numeric(scad_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## MCP
  expt = try(expr = {
    r_mcp14[rpts,] = as.numeric(mcp_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## ISIS
  expt = try(expr = {
    r_isis14[rpts,] = as.numeric(isis_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ##SIS
  expt = try(expr = {
    r_sis14[rpts,] = as.numeric(sis_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  expt = try(expr = {
    r_boruta14[rpts,] = as.numeric(boruta_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## VSURF
  expt = try(expr = {
    r_vsurf14[rpts,] = as.numeric(vsurf_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## RRF
  expt = try(expr = {
    r_rrf14[rpts,] = as.numeric(rrf_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## PIMP
  expt = try(expr = {
    r_pimp14[rpts,] = as.numeric(pimp_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  # ## NTA
  expt = try(expr = {
    r_nta14[rpts,] =  as.numeric(nta_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  ##sivs
  # Stable Iterative Variable Selection (SIVS)
  expt = try(expr = {
    r_sivs14[rpts,] = as.numeric(sivs_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
  },silent = T)
  #
  
  
  ## Stop time
  
  endt <- Sys.time()
  time_taken <- difftime(endt,begin, min)
  print(rpts)
  print(time_taken)
  
  
  
  
}

r_lasso14 = data.frame(r_lasso14)
r_enet14 = data.frame(r_enet14)
r_alasso14 = data.frame( r_alasso14)
r_sivs14 = data.frame(r_sivs14)
r_sparse14 = data.frame(r_sparse14)
r_best14 = data.frame(r_best14)
r_scad14 = data.frame(r_scad14)
r_mcp14 = data.frame(r_mcp14)
r_sis14 = data.frame(r_sis14)
r_isis14 = data.frame(r_isis14)
r_boruta14 = data.frame(r_boruta14)
r_vsurf14 = data.frame(r_vsurf14)
r_rrf14 = data.frame(r_rrf14)
r_pimp14 = data.frame(r_pimp14)
r_nta14 =  data.frame(r_nta14)
colnames(r_lasso14) = colnames(r_enet14) = colnames(r_alasso14) =
  colnames(r_sivs14) =colnames(r_sparse14) =
  colnames(r_best14) =colnames(r_scad14) =colnames(r_mcp14) =
  colnames(r_sis14) =colnames(r_isis14) =colnames(r_boruta14) =
  colnames(r_vsurf14) =colnames(r_rrf14) =colnames(r_pimp14) = colnames(r_nta14) = 
  c("selected","Imp%","Unimp%","MSE","MAE","Accuracy","Precision","Recall","Time","Brier")

#    sim = sim + 1


res_summary14 = data.frame(cbind(colMeans(r_lasso14,na.rm = T),colMeans(r_alasso14,na.rm = T),
                               colMeans(r_sparse14,na.rm = T),colMeans(r_enet14,na.rm = T),
                               colMeans(r_best14,na.rm = T),colMeans(r_sivs14,na.rm = T),
                               colMeans(r_scad14,na.rm = T),colMeans(r_mcp14,na.rm = T),
                               colMeans(r_sis14,na.rm = T),colMeans(r_isis14,na.rm = T),
                               colMeans(r_boruta14,na.rm = T),colMeans(r_vsurf14,na.rm = T),
                               colMeans(r_rrf14,na.rm = T),colMeans(r_pimp14,na.rm = T),
                               colMeans(r_nta14,na.rm = T)))
colnames(res_summary14) = c("LASSO","ALASSO","SparseStep","ElasticNet","Best Subset","SCAD","MCP","SIS","ISIS","SIVS","Boruta","VSURF", "RRF",
                            "PIMP","NTA")
res_summary14 %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()















