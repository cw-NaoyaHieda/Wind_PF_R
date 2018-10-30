
rm(list=ls())
targetPackages <- c("circular","CircStats","grid","reshape2","ggplot2","scales","tidyverse","snow","gridExtra") 
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
# パラメータは真値
phi1 = 0.97;
gam = 3; #constants in wind speed
mu_g = 0.0 ; # location in wind direction for transition
mu_f = 0.0; # location in wind direction for marginal
rho_f = 0.1; # consentration in wind direction for marginal
V = 20;
mu_rho = 0.5;
sig_rho=1
state <- 1 #状態変数の答えが存在するデータなら1



#観測変数の取り出し

y = sample$y;
v = sample$v;
if(state == 1){
  alp = sample$alp;
  rho_ans = sample$r;
}




par1 = c(phi1, gam, mu_g, mu_f, rho_f, V, mu_rho, sig_rho);

#　並列化フィルタリング
X <- particlefilter(par1, y, v, 100)


pfOut1 <- X$pfOut1
rho1 <- X$rho1
wt <- X$wt

# 平滑化
smwt<-particlesmoother(phi, pfOut1, wt)


# 色々Plot thetaなし

# v
mean_v = rowSums(par1[2]*exp(pfOut1/2)*wt)
sm_v = rowSums(par1[2]*exp(pfOut1/2)*smwt)
check_smoothing <- data.frame(dt=1:length(mean_v),mean_v,sm_v,answer = v[-1])

p2 <- ggplot(check_smoothing %>% gather(id,v,-dt), aes(x=dt,y=v,color=id)) +
  geom_line()
print(p2)


# alpha

mean_alpha = rowSums(wt * pfOut1)
sm_alpha = rowSums(smwt * pfOut1)
if(state == 1){
  check_smoothing <- data.frame(dt=1:length(mean_alpha),mean_alpha,sm_alpha,answer = alp[-length(mean_alpha)])
}else{
  check_smoothing <- data.frame(dt=1:length(mean_alpha),mean_alpha,sm_alpha) 
}

p3 <- ggplot(check_smoothing %>% gather(id,alpha,-dt), aes(x=dt,y=alpha,color=id)) +
  geom_line()
print(p3)

mean_rho = rowSums(wt * rho1)
sm_rho = rowSums(smwt * rho1)
if(state == 1){
  check_smoothing <- data.frame(dt=1:length(mean_rho),mean_rho,sm_rho,answer = rho_ans[-length(rho_ans)])
}else{
  check_smoothing <- data.frame(dt=1:length(mean_rho),mean_rho,sm_rho) 
}

# rho

p4 <- ggplot(check_smoothing %>% gather(id,rho,-dt), aes(x=dt,y=rho,color=id)) +
  geom_line()
print(p4)

print(grid.arrange(p2, p3, p4,
                   ncol = 1))

pw_weight <- pairwise_weight(phi1, pfOut1, wt, smwt)



result <- Q_calc_Estep(par1, pfOut1, rho1, pw_weight, smwt, y, v)

print(result)


