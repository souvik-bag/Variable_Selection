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
#       cor_matrix[i,j] = 0.7^(abs(i-j))
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









n = c(100,100,1000)
m = c(10,100,5000)
k = c(3,5,20)
setups = data.frame(n,m,k)
rptns = 100
rpts = 1

i = 3

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


  
  r_lasso = matrix(nrow = rptns,ncol = 10)
  r_enet = matrix(nrow = rptns,ncol = 10)
  r_alasso = matrix(nrow = rptns, ncol=10)
  r_sivs = matrix(nrow = rptns,ncol = 10)
  r_sparse = matrix(nrow = rptns,ncol = 10)
  r_best = matrix(nrow = rptns,ncol = 10)
  r_scad = matrix(nrow = rptns,ncol = 10)
  r_mcp = matrix(nrow = rptns,ncol = 10)
  r_sis = matrix(nrow = rptns,ncol = 10)
  r_isis = matrix(nrow = rptns,ncol = 10)
  r_boruta = matrix(nrow = rptns,ncol = 10)
  r_vsurf = matrix(nrow = rptns,ncol = 10)
  r_rrf = matrix(nrow = rptns,ncol = 10)
  r_pimp =  matrix(nrow = rptns,ncol = 10)
  r_nta =  matrix(nrow = rptns,ncol = 10)
  
  
  
  for (rpts in 1:rptns)
  {
    eps <- rnorm(setups[i,]$n * setups[i,]$m)
    Meps <- matrix(eps, ncol =setups[i,]$n , byrow = TRUE)    
    Meps <- kk %*% Meps
    X <- t(Meps)
    # X= matrix(gen_cor_h(setups[i,]$n,setups[i,]$m,0.5),nrow = setups[i,]$n,ncol = setups[i,]$m)
    
    
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
    ## Lasso
    r_lasso[rpts,] = as.numeric(lasso_simulation(train.data = train.data, test.data =  test.data, beta = beta ))
    ## Elastic net
    r_enet [rpts,] = as.numeric(elastic_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## Awdaptive lasso
    r_alasso[rpts,] = as.numeric( alasso_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## Stable Iterative Variable Selection (SIVS)
    r_sivs[rpts,] = as.numeric(sivs_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## L0L2
    r_best[rpts,] = as.numeric(best_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## Sparse Step
    r_sparse[rpts,] = as.numeric(sparse_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## SCAD
    r_scad[rpts,] = as.numeric(scad_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## MCP
    r_mcp[rpts,] = as.numeric(mcp_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## ISIS
    r_isis[rpts,] = as.numeric(isis_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ##SIS
    r_sis[rpts,] = as.numeric(sis_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## Boruta
    r_boruta[rpts,] = as.numeric(boruta_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## VSURF
    r_vsurf[rpts,] = as.numeric(vsurf_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## RRF
    r_rrf[rpts,] = as.numeric(rrf_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## PIMP
    r_pimp[rpts,] = as.numeric(pimp_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    ## NTA
    r_nta[rpts,] =  as.numeric(nta_simulation(train.data = train.data, test.data =  test.data, beta = beta))
    
    
    
    ## Stop time
    
    endt <- Sys.time()
    time_taken <- difftime(endt,begin, min)
    print(rpts)
    print(time_taken)
    
    
    
  }
  
  
  
  
  r_lasso = data.frame(r_lasso)
  r_enet = data.frame(r_enet)
  r_alasso = data.frame( r_alasso)
  r_sivs = data.frame(r_sivs)
  r_sparse = data.frame(r_sparse)
  r_best = data.frame(r_best)
  r_scad = data.frame(r_scad)
  r_mcp = data.frame(r_mcp)
  r_sis = data.frame(r_sis)
  r_isis = data.frame(r_isis)
  r_boruta = data.frame(r_boruta)
  r_vsurf = data.frame(r_vsurf)
  r_rrf = data.frame(r_rrf)
  r_pimp = data.frame(r_pimp)
  r_nta =  data.frame(r_nta)
  colnames(r_lasso) = colnames(r_enet) = colnames(r_alasso) =
    colnames(r_sivs) =colnames(r_sparse) =
    colnames(r_best) =colnames(r_scad) =colnames(r_mcp) =
    colnames(r_sis) =colnames(r_isis) =colnames(r_boruta) =
    colnames(r_vsurf) =colnames(r_rrf) =colnames(r_pimp) = colnames(r_nta) = 
    c("selected","Imp%","unimp%","mse","mae","misclas","pre","reca","time","brier")
  
  #    sim = sim + 1


res_summary = data.frame(cbind(colMeans(r_lasso,na.rm = T),colMeans(r_alasso,na.rm = T),
                                colMeans(r_sparse,na.rm = T),colMeans(r_enet,na.rm = T),
                                colMeans(r_best,na.rm = T),colMeans(r_sivs,na.rm = T),
                                colMeans(r_scad,na.rm = T),colMeans(r_mcp,na.rm = T),
                                colMeans(r_sis,na.rm = T),colMeans(r_isis,na.rm = T),
                                colMeans(r_boruta,na.rm = T),colMeans(r_vsurf,na.rm = T),
                                colMeans(r_rrf,na.rm = T),colMeans(r_pimp,na.rm = T),
                                colMeans(r_nta,na.rm = T)))
colnames(res_summary) = c("lasso","alasso","sparse","enet","best","sivs","scad","mcp","sis","isis","boruta","vsurf", "rrf",
                           "pimp","nta")
res_summary %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()
