
rm(list=ls())
targetPackages <- c("circular","CircStats","grid","reshape2","ggplot2","scales","tidyverse","snow") 
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
source("functions/WSD_functions.R")
source("functions/function.R")
set.seed(1000)
#時系列数
N = 800;

#塩浜先生がお作りしたものを使用しています
sample = read.csv("sample_1000.csv");
sample = sample[1:800,];
#観測変数の取り出し(観測変数としては使わないものもある)

y = sample$y;
v = sample$v;

opt_params <- matrix(NA,ncol=11,nrow=100)
#パラメータ設定 微妙にずらします

i = 0
hosts <- rep('localhost',12)
scl <- makeCluster(hosts, "SOCK")
for(j in 1:5){
  
  phi1 = min(0.97+rnorm(1,sd=0.01),0.99);
  gam = max(3+rnorm(1),1); #constants in wind speed
  mu_g = 0.0 + rnorm(1,0.5); # location in wind direction for transition
  mu_f = 0.0+ rnorm(1,0.5); # location in wind direction for marginal
  rho_f = max(0.1+ rnorm(1,0.05),0.01); # consentration in wind direction for marginal
  V = 20+rnorm(1,sd=2);
  mu_rho = 0.5+rnorm(1,sd=0.2);19931218ab
  sig_rho=max(1+rnorm(1,sd=0.3),0.01);
  par1 = c(phi1, gam, mu_g, mu_f, rho_f, V, mu_rho, sig_rho);
  
  #初期パラメータの設定
  
  opt_params[i+1,]= c(j,1,0,par1)
  i <- i+1
  
  k=2
  while(k==2 || (check/800 > 0.0001 && k<=11) ){
    X <- particlefilter(par1, y, v, 100)
    pfOut1 <- X$pfOut1
    rho1 <- X$rho1
    wt <- X$wt
    phi <- par1[1]
    smwt<-particlesmoother(phi, pfOut1, wt)
    
    #mean_alpha = rowSums(wt * pfOut1)
    #sm_alpha = rowSums(smwt * pfOut1)
    #check_smoothing <- data.frame(dt=1:799,mean_alpha[],sm_alpha[])
    
    pw_weight <- pairwise_weight(phi1, pfOut1, wt, smwt)
    
    par1[1] = sig_env(par1[1])
    par1[2] = log(par1[2])
    par1[5] = sig_env(par1[5])
    par1[6] = log(par1[6])
    par1[8] = log(par1[8])
    
    
   
    
    
    
    optim_fun <- function(par1){
      Q_calc_para(par1, pfOut1, rho1, pw_weight, smwt, y, v)
    }
    result <- optim(par = par1,fn = optim_fun)
    
    
    par1[1] = sig(result$par[1])
    par1[2] = exp(result$par[2])
    par1[5] = sig(result$par[5])
    par1[6] = exp(result$par[6])
    par1[8] = exp(result$par[8])
    opt_params[i+1,] <- c(j,k,result$value,par1)
    i <- i+1
    check <- abs(opt_params[i,3] - opt_params[i-1,3])
    print(i)
    print(j)
    print(k)
    print(check)
    k <- k+1
  }
}

stopCluster(scl)