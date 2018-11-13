rm(list=ls())
Packages <- c("circular","CircStats","reshape2","scales","tidyverse","optimParallel","snow") 
for(package in Packages) library(package, character.only = T)
source("functions/WSD_functions.R")
source("functions/function.R")
set.seed(1000)

phi1 = 0.90
gam = 3 #constants in wind speed
mu_g = 0.0  # location in wind direction for transition
mu_f = 0.0 # location in wind direction for marginal
rho_f = 0.1 # consentration in wind direction for marginal
V = 3
mu_rho = 0.5
sig_rho= 1
par1 = c(phi1, gam, mu_g, mu_f, rho_f, V, mu_rho, sig_rho);

x <-seq(0.001, 8, by=0.01)
y_slide <- (tanh(sig_rho*x+mu_rho)+1)/2
ggplot(data.frame(x,y_slide) %>% gather(id,value,-x),aes(x,value,color=id))+geom_line()

z <- simulate.data(800, par=par1)
y <- z$theta
alp <- z$alpha
v <- z$v
n <- length(y)

ggplot(data.frame(y,v) ,aes(x=v,y=y))+geom_point()

df=data.frame(x=1:n,
              y=z$theta,
              v=z$v,
              rho=z$rho,
              alpha=z$alpha)

df2=melt(df,id.vars=c("x"),
         measure.vars=c("y","v","rho", "alpha"))
levels(df2$variable)[levels(df2$variable)=="y"] <- "Wind Direction"
levels(df2$variable)[levels(df2$variable)=="v"] <- "Wind Speed"

ggplot(df2,aes(x,value))+
  facet_grid(variable~.,scales="free_y")+ylab("")+
  geom_point(data = subset(df2, variable == "Wind Direction"),size=1.2) +
  geom_line(data = subset(df2, variable == "Wind Speed"),size=.1) +
  geom_line(data = subset(df2, variable == "rho"),size=.1)+
  geom_line(data = subset(df2, variable == "alpha"),size=.1)



X_100 <- particlefilter(par1, y, v, 100)
pfOut1_100 <- X_100$pfOut1
rho1_100 <- X_100$rho1
wt_100 <- X_100$wt
setwd("~/Desktop/Wind_PF_R/shiny")

write.csv(df,"df.csv")
write.csv(pfOut1_100 ,'pfOut1_100.csv')
write.csv(rho1_100  ,'rho1_100.csv')
write.csv(wt_100,'wt_100.csv')

X_1000 <- particlefilter(par1, y, v, 1000)
pfOut1_1000 <- X_1000$pfOut1
rho1_1000 <- X_1000$rho1
wt_1000 <- X_1000$wt
write.csv(pfOut1_1000 ,'pfOut1_1000.csv')
write.csv(rho1_1000  ,'rho1_1000.csv')
write.csv(wt_1000,'wt_1000.csv')


X_10000 <- particlefilter(par1, y, v, 10000)
pfOut1_10000 <- X_10000$pfOut1
rho1_10000 <- X_10000$rho1
wt_10000 <- X_10000$wt
write.csv(pfOut1_10000 ,'pfOut1_10000.csv')
write.csv(rho1_10000  ,'rho1_10000.csv')
write.csv(wt_10000,'wt_10000.csv')

mean_v_100 = rowSums(par1[2]*exp(pfOut1_100/2)*wt_100)
mean_v_1000 = rowSums(par1[2]*exp(pfOut1_1000/2)*wt_1000)
mean_v_10000 = rowSums(par1[2]*exp(pfOut1_10000/2)*wt_10000)
check_mean <- data.frame(dt=1:length(mean_v_100),mean_v_100,mean_v_1000, mean_v_10000,answer = v[-length(v)])

p2 <- ggplot(check_mean %>% gather(id,v,-dt), aes(x=dt,y=v,color=id)) +
  geom_line()
print(p2)


mean_alpha_100 = rowSums(wt_100 * pfOut1_100)
mean_alpha_1000 = rowSums(wt_1000 * pfOut1_1000)
mean_alpha_10000 = rowSums(wt_10000 * pfOut1_10000)

check_mean <- data.frame(dt=1:length(mean_alpha_100),mean_alpha_100,
                              mean_alpha_1000,
                              mean_alpha_10000,
                              answer = alp[-length(mean_alpha_100)])

p3 <- ggplot(check_mean %>% gather(id,alpha,-dt), aes(x=dt,y=alpha,color=id)) +
  geom_line()
print(p3)


plot(pfOut1_100[200,],wt_100[200,])
plot(pfOut1_1000[200,],wt_1000[200,])
plot(pfOut1_10000[200,],wt_10000[200,])
tmp <- data.frame(x = round(pfOut1_100[200,],digits = 2), y = wt_100[200,])
ggplot(tmp,aes(x,y))+geom_bar(stat = "identity") +
  geom_vline(xintercept = sum(wt_100[200,] * pfOut1_100[200,]),color='blue')
tmp <- data.frame(x = round(pfOut1_1000[200,],digits = 2), y = wt_1000[200,])
ggplot(tmp,aes(x,y))+geom_bar(stat = "identity")+
  geom_vline(xintercept = sum(wt_1000[200,] * pfOut1_1000[200,]),color='blue')
tmp <- data.frame(x = round(pfOut1_10000[200,],digits = 2), y = wt_10000[200,])
ggplot(tmp,aes(x,y))+geom_bar(stat = "identity")
