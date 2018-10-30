
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
#観測変数の取り出し

y = sample$y;
v = sample$v;



phi1 = 0.97;
gam = 3; #constants in wind speed
mu_g = 0.0 ; # location in wind direction for transition
mu_f = 0.0; # location in wind direction for marginal
rho_f = 0.1; # consentration in wind direction for marginal
V = 20;
mu_rho = 0.5;
sig_rho=1
par1 = c(phi1, gam, mu_g, mu_f, rho_f, V, mu_rho, sig_rho);
  
#初期パラメータの設定
  
opt_params[i+1,]= c(j,1,0,par1)

X <- particlefilter(par1, y, v, 100)

pfOut1 <- X$pfOut1
rho1 <- X$rho1
wt <- X$wt

smwt<-particlesmoother(phi, pfOut1, wt)
    
    #mean_alpha = rowSums(wt * pfOut1)
    #sm_alpha = rowSums(smwt * pfOut1)
    #check_smoothing <- data.frame(dt=1:799,mean_alpha[],sm_alpha[])
    
pw_weight <- pairwise_weight(phi1, pfOut1, wt, smwt)
    

result <- Q_calc(par1, pfOut1, rho1, pw_weight, smwt, y, v)

result


