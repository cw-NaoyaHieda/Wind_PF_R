###
### wrapped Cauchy distribution??score function
###


score.wrappedcauchy <- function(y, mu, rho)
{
  tmp1 <- - 2 * rho * sin(y-mu)
  tmp2 <- 1+rho^2-2*rho*cos(y-mu)
  return(tmp1/tmp2)
}


#dth <- seq(-pi, pi, length=1000)
#out <- sapply(dth, score.wrappedcauchy, mu=pi/4, rho=0.4)
#plot(dth, out, type="l")
#lines(dth,dwrpcauchy(dth, mu=pi/4, rho=0.4),col=2)


###
### psswrapedcauchy: sswc
###

psswrappedcauchy <-function (y, mu , rho, lambda )
{
  ##lambda=0 wrapped cauchy
  ## -pi < y < pi
  tmp1 <- (1+rho)/(1-rho)*tan( (y-mu)/2 )
  tmp2 <- lambda*(1-rho^2)/(4*pi*rho)*log( (1+rho^2-2*rho*cos(y-mu))/(1+rho)^2 )
  tmp <-unname(0.5+ atan(tmp1)/pi+tmp2)
  return(tmp)
}



###
### dsswrapedcauchy: sswc
###

dsswrpcauchy<- function(theta, mu, rho, lambda){
  (1 - rho^2)/((2 * pi) * (1 + rho^2 - 2 * rho * cos(theta - 
                                                       mu)))*(1+lambda*sin(theta-mu))
}


###
### d_conditional_WJ:
###

d_conditional_WJ <- function(x, theta0, mu_g, rho_g,  mu_f, rho_f, q=1){
  tmp1 <- psswrappedcauchy(x, mu=mu_f, rho=rho_f,lambda=0)
  tmp2 <- psswrappedcauchy(theta0, mu=mu_f, rho=rho_f,lambda=0) 
  out1 <- dsswrpcauchy( 2*pi*(tmp1-q*tmp2), mu =mu_g, rho=rho_g,lambda=0)
  out <- 2*pi*out1*dsswrpcauchy(x, mu =mu_f, rho=rho_f,lambda=0)
  return(out) 
}


###
### d_conditional_WJ
###

#library(ggplot2)
#f <- ggplot(data.frame(x = c(-pi, pi)), aes(x))
#f + stat_function(fun = d_conditional_WJ, colour = "black",n=400,
#                  args=list(theta0=0,mu_g=0.01, rho_g=0.8, mu_f=0.1, rho_f=0.4)) +
#    stat_function(fun = d_conditional_WJ, colour = "red", n=400,
#                  args=list(theta0=0,mu_g=-0.01, rho_g=0.4, mu_f=0.2, rho_f=0.6)) +
#    stat_function(fun = d_conditional_WJ, colour = "blue", n=400,
#                  args=list(theta0=0,mu_g=0.0, rho_g=0.6, mu_f=-0.3, rho_f=0.8)) 
###
###
###  d_conditional_WJ
###

r_conditional_WJ <- function(n, theta0, mu_g, rho_g,  mu_f, rho_f, q=1){
  dx <- seq(-pi,pi,by=0.01)
  M <- max(sapply(dx, d_conditional_WJ,
                  theta0=theta0,mu_g=mu_g, rho_g=rho_g, mu_f=mu_f, rho_f=rho_f, q=1))
  i <- 1
  result <- c(1:n)
  while (i <= n) {
    x <- runif(1, -pi,  pi)
    y <- M*runif(1,0,1)
    f <- d_conditional_WJ(x,theta0, mu_g, rho_g,  mu_f, rho_f)
    if (y <= f) {
      result[i] <- x
      i <- i + 1
    }
  }
  result
}


### pforeach
### [pforeach](http://d.hatena.ne.jp/hoxo_m/20141222/p1)
###


#library(pforeach)
r_conditional_WJ2 <- function(n, theta0, mu_g, rho_g,  mu_f, rho_f, q=1){
  dx <- seq(-pi,pi,by=0.01)
  M <- max(sapply(dx, d_conditional_WJ,
                  theta0=theta0,mu_g=mu_g, rho_g=rho_g, mu_f=mu_f, rho_f=rho_f, q=1))
  i <- 1
  result <- c(1:n)
  result <- pforeach(i = 1 : n)({
    reject=TRUE
    while( reject){
      x <- runif(1, -pi,  pi)
      y <- M*runif(1,0,1)  
      f <- d_conditional_WJ(x,theta0, mu_g, rho_g,  mu_f, rho_f)  
      if (y <= f) { 
        out =x
        reject=FALSE}}
    out
  })
}

r_conditional_WJ_mean <- function(n, theta0, mu_g, rho_g,  mu_f, rho_f, q=1){
  dx <- seq(-pi,pi,by=0.01)
  M <- max(sapply(dx, d_conditional_WJ,
                  theta0=theta0,mu_g=mu_g, rho_g=rho_g, mu_f=mu_f, rho_f=rho_f, q=1))
  i <- 1
  result <- c(1:n)
  while (i <= n) {
    x <- runif(1, -pi,  pi)
    y <- M*runif(1,0,1)
    f <- d_conditional_WJ(x,theta0, mu_g, rho_g,  mu_f, rho_f)
    if (y <= f) {
      result[i] <- x
      i <- i + 1
    }
  }
  r_s <- sin(result)
  r_c <- cos(result)
  
  return(circular.mean(mean(r_s),mean(r_c)))
}


##

##
#t<-proc.time()
#x.sim1<-r_conditional_WJ(n=2000, 
#          theta0=0, mu_g=0.01, rho_g=.8,  mu_f=0.1, rho_f=0.4, q=1)
#proc.time()-t
##
#t<-proc.time()
#x.sim2<-r_conditional_WJ2(n=2000, 
#          theta0=0, mu_g=0.01, rho_g=.8,  mu_f=0.1, rho_f=0.4, q=1)
#proc.time()-t




#x.sim1<-r_conditional_WJ(n=1000, 
#                        theta0=0, mu_g=0.01, rho_g=.8,  mu_f=0.1, rho_f=0.4, q=1)
#df <- data.frame(x=x.sim1)
#ggplot(df, aes(x))+
#        geom_histogram(binwidth=.05, colour="black", 
#                          aes(y=..density.., fill=..count..))+
#       scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
#       xlim(-pi,pi)+
#       stat_function(fun = d_conditional_WJ, colour = "red",n=400,
#                  args=list(theta0=0,mu_g=0.01, rho_g=0.8, mu_f=0.1, rho_f=0.4))
#x.sim1<-r_conditional_WJ(n=1000, 
#                        theta0=pi/2, mu_g=0.01, rho_g=.2,  mu_f=0.3, rho_f=0.8, q=1)
#df <- data.frame(x=x.sim1)
#ggplot(df, aes(x))+
#        geom_histogram(binwidth=.05, colour="black", 
#                          aes(y=..density.., fill=..count..))+
#       scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
#       xlim(-pi,pi)+
#       stat_function(fun = d_conditional_WJ, colour = "red",n=400,
#                  args=list(theta0=pi/2,mu_g=0.01, rho_g=0.2, mu_f=0.3, rho_f=0.8))



###
###
###

simulate.data <- function(n, par){
  # n:sample size
  # par: parameter
  # ini: means and variances of initial state 
  phi1 <- par[1] # AR in state of wind speed 
  gam  <- par[2] # constants in log wind speed 
  mu_g <- par[3] # location in wind direction for transition  
  mu_f <- par[4] # location in wind direction for marginal
  rho_f <-par[5] # consentration in wind direction for marginal
  V <- par[6]
  mu_rho <- par[7]
  sig_rho <- par[8]
  ## initial state
  theta <- rwrpcauchy(1, location=mu_f, rho=rho_f)
  theta <- ifelse(theta>pi, theta-2*pi, theta)
  
  alpha <- 0
  v  <- gam * rgamma(1, shape=V, rate=V)
  
  for(i in 2:n){ 
    
    alpha[i] <- phi1*alpha[i-1] +sqrt(1-phi1^2)*rnorm(1)
    theta[i] <- r_conditional_WJ(n=1, 
                                 theta0=theta[i-1], mu_g=mu_g,   
                                 rho_g=0.95*(tanh(sig_rho*alpha[i-1]+mu_rho)+1)/2,
                                 mu_f=mu_f, rho_f=rho_f, q=1)
    
    theta[i] <- ifelse(theta[i] > pi, theta[i]-2*pi, theta[i])
    theta[i] <- ifelse(theta[i] < -pi, theta[i]+2*pi, theta[i])
    #print(0.95*(tanh(sig_rho*alpha[i-1]+mu_rho)+1)/2)
    
    v[i] <-   gam*exp(alpha[i]/2)*rgamma(1, shape=V, rate=V)
  }
  rho =0.95*(tanh(sig_rho*alpha+mu_rho)+1)/2
  return(list(alpha=as.numeric(alpha),
              theta=as.numeric(theta),
              v=as.numeric(v), 
              rho=as.numeric(rho)))
}




## Miscellaneous  function

circular.mean <- function(m.sin, m.cos){
  if(m.sin >0 &  m.cos> 0) return(atan(m.sin/m.cos))
  else if(m.sin < 0 &  m.cos >  0) return(atan(m.sin/m.cos))
  else if(m.sin < 0 &  m.cos <  0) return(atan(m.sin/m.cos)-pi)
  else if(  m.sin >0 & m.cos <  0) return(atan(m.sin/m.cos)+pi)
}


Resample1 <- function(data, weight, NofSample){
  re_ind <- runif(NofSample)
  cmwt <- cumsum(as.double(weight));
  st <- sapply(re_ind, function(x) sum(x>cmwt[-length(cmwt)]))
  newdata <- data[ (st+1) ]
  return(newdata)
}


## Filtering function


Resample2 <- function(data1, data2, weight, NofSample){
  re_ind <- runif(NofSample)
  cmwt <- cumsum(weight)/sum(weight);
  st <- sapply(re_ind, function(x) sum(x>cmwt[-length(cmwt)]))
  newdata1 <- data1[ (st+1) ]
  newdata2 <- data2[ (st+1) ]
  return(list(newdata1, newdata2) )
}



Resample3 <- function(data1, data2, weight, NofSample){
  fit1 <- density(data1, window="epanechnikov", weights =weight)
  data1.new <- rnorm(NofSample, sample(data1, size = NofSample, replace = TRUE), fit1$bw)
  fit2 <- density(data2, window="epanechnikov", weights =weight) 
  data2.new <- rnorm(NofSample, sample(data2, size = NofSample, replace = TRUE), fit2$bw)
  return(list(data1=data1.new, data2=data2.new) )
}

Resample_gpu <- function(data, weight, NofSample){
  re_ind <- gpuVector(runif(NofSample))
  cmwt <- cumsum(weight[]);
  st <- sapply(re_ind, function(x) sum(x>cmwt[-length(cmwt)]))
  newdata <- data[][c(st+1) ]
  return(newdata)
}

sig <- function(x){
  y = (tanh(x) + 1) / 2;
  return(y)
}
sig_env <- function(y){
  x = (1/2)*log(y/(1-y));
  return(x)
}

pi_shori <- function (x){
  while(x> pi){
    x = x-2*pi; 
  }
  
  while(x < -pi){
    x = x+2*pi;
    end
  }
  y = x;
  
  return(y)
}


sig_95 <- function(x){
  y = 0.95*(tanh(x) + 1) / 2;
  return(y)
}
sig_env_95 <- function(y){
  x = (1/2)*log(2*y/(1.9-2*y));
  return(x)
}