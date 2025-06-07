
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


## create setup frame
n = c(100,100,1000)
m = c(10,100,5000)
k = c(3,5,20)
setups = data.frame(n,m,k)
rptns = 100
rpts = 1


# Select the value of i

i <- 1 

r_lasso11 = matrix(nrow = rptns,ncol = 10)
r_enet11 = matrix(nrow = rptns,ncol = 10)
r_alasso11 = matrix(nrow = rptns, ncol=10)
r_sivs11 = matrix(nrow = rptns,ncol = 10)
r_sparse11 = matrix(nrow = rptns,ncol = 10)
r_best11 = matrix(nrow = rptns,ncol = 10)
r_scad11 = matrix(nrow = rptns,ncol = 10)
r_mcp11 = matrix(nrow = rptns,ncol = 10)
r_sis11 = matrix(nrow = rptns,ncol = 10)
r_isis11 = matrix(nrow = rptns,ncol = 10)
r_boruta11 = matrix(nrow = rptns,ncol = 10)
r_vsurf11 = matrix(nrow = rptns,ncol = 10)
r_rrf11 = matrix(nrow = rptns,ncol = 10)
r_pimp11 =  matrix(nrow = rptns,ncol = 10)
r_nta11 =  matrix(nrow = rptns,ncol = 10)



for (rpts in 1:rptns)
{
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
  ## Lasso
  r_lasso11[rpts,] = as.numeric(lasso_simulation(train.data = train.data, test.data =  test.data, beta = beta ))
  ## Elastic net
  r_enet11[rpts,] = as.numeric(elastic_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## Awdaptive lasso
  r_alasso11[rpts,] = as.numeric( alasso_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## Stable Iterative Variable Selection (SIVS)
  r_sivs11[rpts,] = as.numeric(sivs_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## L0L2
  r_best11[rpts,] = as.numeric(best_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## Sparse Step
  r_sparse11[rpts,] = as.numeric(sparse_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## SCAD
  r_scad11[rpts,] = as.numeric(scad_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## MCP
  r_mcp11[rpts,] = as.numeric(mcp_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## ISIS
  r_isis11[rpts,] = as.numeric(isis_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ##SIS
  r_sis11[rpts,] = as.numeric(sis_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## Boruta
  r_boruta11[rpts,] = as.numeric(boruta_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## VSURF
  r_vsurf11[rpts,] = as.numeric(vsurf_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## RRF
  r_rrf11[rpts,] = as.numeric(rrf_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## PIMP
  r_pimp11[rpts,] = as.numeric(pimp_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  ## NTA
  r_nta11[rpts,] =  as.numeric(nta_simulation(train.data = train.data, test.data =  test.data, beta = beta))
  
  
  
  ## Stop time
  
  endt <- Sys.time()
  time_taken <- difftime(endt,begin, min)
  print(rpts)
  print(time_taken)
  
  
  
}


r_lasso11 = data.frame(r_lasso11)
r_enet11 = data.frame(r_enet11)
r_alasso11 = data.frame( r_alasso11)
r_sivs11 = data.frame(r_sivs11)
r_sparse11 = data.frame(r_sparse11)
r_best11 = data.frame(r_best11)
r_scad11 = data.frame(r_scad11)
r_mcp11 = data.frame(r_mcp11)
r_sis11 = data.frame(r_sis11)
r_isis11 = data.frame(r_isis11)
r_boruta11 = data.frame(r_boruta11)
r_vsurf11 = data.frame(r_vsurf11)
r_rrf11 = data.frame(r_rrf11)
r_pimp11 = data.frame(r_pimp11)
r_nta11 =  data.frame(r_nta11)
colnames(r_lasso11) = colnames(r_enet11) = colnames(r_alasso11) =
  colnames(r_sivs11) =colnames(r_sparse11) =
  colnames(r_best11) =colnames(r_scad11) =colnames(r_mcp11) =
  colnames(r_sis11) =colnames(r_isis11) =colnames(r_boruta11) =
  colnames(r_vsurf11) =colnames(r_rrf11) =colnames(r_pimp11) = colnames(r_nta11) = 
  c("selected","Imp%","unimp%","mse","mae","misclas","pre","reca","time","brier")

#    sim = sim + 1


res_summary11 = data.frame(cbind(colMeans(r_lasso11,na.rm = T),colMeans(r_alasso11,na.rm = T),
                               colMeans(r_sparse11,na.rm = T),colMeans(r_enet11,na.rm = T),
                               colMeans(r_best11,na.rm = T),colMeans(r_sivs11,na.rm = T),
                               colMeans(r_scad11,na.rm = T),colMeans(r_mcp11,na.rm = T),
                               colMeans(r_sis11,na.rm = T),colMeans(r_isis11,na.rm = T),
                               colMeans(r_boruta11,na.rm = T),colMeans(r_vsurf11,na.rm = T),
                               colMeans(r_rrf11,na.rm = T),colMeans(r_pimp11,na.rm = T),
                               colMeans(r_nta11,na.rm = T)))
colnames(res_summary11) = c("lasso","alasso","sparse","enet","best","sivs","scad","mcp","sis","isis","boruta","vsurf", "rrf",
                          "pimp","nta")
res_summary11 %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()















