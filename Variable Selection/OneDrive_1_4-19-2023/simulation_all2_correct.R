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

## correlated data generation
# Generate hierarchical Correlation 

# cor_matrix = matrix(nrow = m, ncol = m)
# diag(cor_matrix) = 1
# 
# i =1
# j= 1
# for (i in 1:m)
#   {
#   
#  for (j in 1:m)
#    {
# 
#     cor_matrix[i,j] = 0.7^(abs(i-j))
#     j <- j+1
#  }
#   i <- i+1
#   }
# 
# cor_matrix

# gen_cor_h <- function(n,m,cor){
#   
#   cor_matrix <- matrix(nrow = m, ncol = m)
#   diag(cor_matrix) <- 1
#   
#   # generate hierarchical correlation
#   i =1
#   j= 1
#   for (i in 1:m)
#   {
#     
#     for (j in 1:m)
#     {
#       
#       cor_matrix[i,j] = 0.5^(abs(i-j))
#       j <- j+1
#     }
#     i <- i+1
#   }
#   
#   
#   SigmaEV <- eigen(cor_matrix)
#   eps <- rnorm(n * ncol(SigmaEV$vectors))
#   Meps <- matrix(eps, ncol = n, byrow = TRUE)    
#   Meps <- SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)) %*% Meps
#   Meps <- t(Meps)
#   return(Meps)
# }


# q = gen_cor_h(1000,10,0.7)
# 
# cor(q)









n = c(100,100,100,100)
m = c(10,100,1000,1000)
k = c(3,5,10,200)
# cor = c(0,0.2,NA)
setups = data.frame(n,m,k)
# setups[1,]$n
# sim = 1
rptns = 100
rpts = 1
#case 1 , n=100, m = 10

i = 4
## Multivariate normal data generation with correlation

cor_matrix <- matrix(nrow = setups[i,]$m, ncol = setups[i,]$m)
diag(cor_matrix) <- 1

# generate hierarchical correlation
r =1
c = 1
for (r in 1:setups[i,]$m)
{
  
  for (c in 1:setups[i,]$m)
  {
    
    cor_matrix[r,c] = 0.5^(abs(r-c))
    c <- c+1
  }
  r <- r+1
}


SigmaEV <- eigen(cor_matrix)

kk <- SigmaEV$vectors %*% diag(sqrt(SigmaEV$values))



  # 
  r_lasso24 = matrix(nrow = rptns,ncol = 10)
  r_enet24 = matrix(nrow = rptns,ncol = 10)
  r_alasso24 = matrix(nrow = rptns, ncol=10)
  r_sivs24 = matrix(nrow = rptns,ncol = 10)
  r_sparse24 = matrix(nrow = rptns,ncol = 10)
  r_best24 = matrix(nrow = rptns,ncol = 10)
  r_scad24 = matrix(nrow = rptns,ncol = 10)
  r_mcp24 = matrix(nrow = rptns,ncol = 10)
  r_sis24 = matrix(nrow = rptns,ncol = 10)
  r_isis24 = matrix(nrow = rptns,ncol = 10)
  r_boruta24 = matrix(nrow = rptns,ncol = 10)
  r_vsurf24 = matrix(nrow = rptns,ncol = 10)
  r_rrf24 = matrix(nrow = rptns,ncol = 10)
  r_pimp24 =  matrix(nrow = rptns,ncol = 10)
  r_nta24 =  matrix(nrow = rptns,ncol = 10)
  
  
  
  for(rpts in 1:rptns)
  {
    set.seed(rpts)
    eps <- rnorm(setups[i,]$n * setups[i,]$m)
    Meps <- matrix(eps, ncol =setups[i,]$n , byrow = TRUE)    
    Meps <- kk %*% Meps
    X <- t(Meps)
    # X= matrix(gen_cor_h(setups[i,]$n,setups[i,]$m,0.5),nrow = setups[i,]$n,ncol = setups[i,]$m)
    # 
    
    beta <- numeric(setups[i,]$m)
    beta[sample(setups[i,]$m,setups[i,]$k)] = rnorm(setups[i,]$k,0,3)
    beta <- ifelse((beta>0 & beta <0.5),beta + 0.5, ifelse((beta <0 & beta > -0.5), beta-0.5,beta))
    #beta[1:setups[i,]$k] = rnorm(setups[i,]$k,0,3)
    #beta <- ifelse((beta>=0 & beta <0.5),beta + 0.5, ifelse((beta <0 & beta > -0.5), beta-0.5,beta))
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
      r_lasso24[rpts,] = as.numeric(lasso_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta ))
    },silent = T)
    ## Elastic net
    expt = try(expr = {
      r_enet24 [rpts,] = as.numeric(elastic_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    ## Awdaptive lasso
    expt = try(expr = {
      r_alasso24[rpts,] = as.numeric( alasso_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)

    # ## L0L2
    expt = try(expr = {
      r_best24[rpts,] = as.numeric(best_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## Sparse Step
    expt = try(expr = {
      r_sparse24[rpts,] = as.numeric(sparse_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## SCAD
    expt = try(expr = {
      r_scad24[rpts,] = as.numeric(scad_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## MCP
    expt = try(expr = {
      r_mcp24[rpts,] = as.numeric(mcp_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## ISIS
    expt = try(expr = {
      r_isis24[rpts,] = as.numeric(isis_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ##SIS
    expt = try(expr = {
      r_sis24[rpts,] = as.numeric(sis_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    expt = try(expr = {
      r_boruta24[rpts,] = as.numeric(boruta_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## VSURF
    expt = try(expr = {
      r_vsurf24[rpts,] = as.numeric(vsurf_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## RRF
    expt = try(expr = {
      r_rrf24[rpts,] = as.numeric(rrf_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## PIMP
    expt = try(expr = {
      r_pimp24[rpts,] = as.numeric(pimp_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    # ## NTA
    expt = try(expr = {
      r_nta24[rpts,] =  as.numeric(nta_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    ##sivs
    # Stable Iterative Variable Selection (SIVS)
    expt = try(expr = {
      r_sivs24[rpts,] = as.numeric(sivs_simulation(X.train = X.train,Y.train =Y.train,X.test = X.test,Y.test = Y.test, beta = beta))
    },silent = T)
    #

    
    ## Stop time
    
    endt <- Sys.time()
    time_taken <- difftime(endt,begin, min)
    print(rpts)
    print(time_taken)
          
    
    
    
  }
  r_lasso24 = data.frame(r_lasso24)
  r_enet24 = data.frame(r_enet24)
  r_alasso24 = data.frame( r_alasso24)
  r_sivs24 = data.frame(r_sivs24)
  r_sparse24 = data.frame(r_sparse24)
  r_best24 = data.frame(r_best24)
  r_scad24 = data.frame(r_scad24)
  r_mcp24 = data.frame(r_mcp24)
  r_sis24 = data.frame(r_sis24)
  r_isis24 = data.frame(r_isis24)
  r_boruta24 = data.frame(r_boruta24)
  r_vsurf24 = data.frame(r_vsurf24)
  r_rrf24 = data.frame(r_rrf24)
  r_pimp24 = data.frame(r_pimp24)
  r_nta24 =  data.frame(r_nta24)
  colnames(r_lasso24) = colnames(r_enet24) = colnames(r_alasso24) =
    colnames(r_sivs24) =colnames(r_sparse24) =
    colnames(r_best24) =colnames(r_scad24) =colnames(r_mcp24) =
    colnames(r_sis24) =colnames(r_isis24) =colnames(r_boruta24) =
    colnames(r_vsurf24) =colnames(r_rrf24) =colnames(r_pimp24) = colnames(r_nta24) = 
    c("selected","Imp%","Unimp%","MSE","MAE","Accuracy","Precision","Recall","Time","Brier")
  
  #    sim = sim + 1


res_summary24 = data.frame(cbind(colMeans(r_lasso24,na.rm = T),colMeans(r_alasso24,na.rm = T),
                                colMeans(r_sparse24,na.rm = T),colMeans(r_enet24,na.rm = T),
                                colMeans(r_best24,na.rm = T),
                                colMeans(r_scad24,na.rm = T),colMeans(r_mcp24,na.rm = T),
                                colMeans(r_sis24,na.rm = T),colMeans(r_isis24,na.rm = T),colMeans(r_sivs24,na.rm = T),
                                colMeans(r_boruta24,na.rm = T),colMeans(r_vsurf24,na.rm = T),
                                colMeans(r_rrf24,na.rm = T),colMeans(r_pimp24,na.rm = T),
                                colMeans(r_nta24 ,na.rm = T)))
colnames(res_summary24) = c("LASSO","ALASSO","SparseStep","ElasticNet","Best Subset","SCAD","MCP","SIS","ISIS","SIVS","Boruta","VSURF", "RRF",
                           "PIMP","NTA")
res_summary24 %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()





