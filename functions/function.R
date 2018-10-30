particlefilter <- function(par, y, v, nParticle){
  
  #パラメータの取得
  print(par)
  phi1 <- par[1]; 
  gam  <- par[2]; 
  mu_g <- par[3];
  mu_f <- par[4]; 
  rho_f <- par[5];
  V <- par[6];
  mu_rho <- par[7];
  sig_rho <- par[8];
  
  dT <- length(y);
  
  
  pfOut1 <- (matrix(nrow=(dT-1), ncol=nParticle));
  
  wt <- (matrix(nrow=(dT-1), ncol=nParticle));
  rho1 <- (matrix(nrow=(dT-1), ncol=nParticle));
  
  a0 = rnorm(nParticle);
  t0 = rwrpcauchy(nParticle, mu_f, rho_f);
  t0 = sapply(t0, pi_shori);
  
  
  pfOut1[1, ] <- a0;
  rho1[1,] <-  0.95*(tanh(sig_rho*a0 + mu_rho )+1)/2;
  wt[1,] <- rep(1 / nParticle,nParticle);
  
  
  
  N_eff = rep(0,nParticle)
  nEff = nParticle/10
  
  
  for(dt in 2:(dT-1)){
    
    pfOut1[dt,] <- phi1 * pfOut1[dt - 1, ] + rnorm(nParticle,sd=sqrt(1-phi1^2))
    rho1[dt,] <- 0.95 * ( tanh( sig_rho * pfOut1[dt,] + mu_rho)+1) / 2
    if(sum(rho1[dt,] == 0)>0){
      rho1[dt,which(rho1[dt,] == 0)] <-  5.2736e-17
    }
    
    tmp1 = d_conditional_WJ(y[dt+1], y[dt], mu_g, rho1[dt,], mu_f, rho_f, 1)
    tmp2 = dgamma(v[dt]/(gam*exp(pfOut1[dt,]/2)) , shape = V, rate = V)/(gam*exp(pfOut1[dt,]/2))
    
    wt[dt,] = (tmp1/sum(tmp1)) * (tmp2/sum(tmp2)) * wt[dt-1,]
    wt_tmp <- wt[dt,] / sum(wt[dt,])
    wt[dt,] <- wt_tmp
    
    N_eff[dt] = 1 / sum(wt[dt,]^2);
    if(is.nan(sum(wt[dt,]))){
      browser()
      print(dt)
    }
    if(N_eff[dt] < nEff){
      pfOut1[dt,] <- Resample1(pfOut1[dt,], wt[dt,], nParticle);
      rho1[dt,] <- Resample1(rho1[dt,], wt[dt,], nParticle);
      wt[dt,] <- rep(1 / nParticle,nParticle);
    }
    if(sum(is.na(wt[dt,]) > 0)){
      browser()
    }
  }
  
  return(list(pfOut1 = pfOut1,rho1 = rho1,wt = wt) )
  
  
}

particlesmoother <- function(phi, pfOut1, wt){
  size = dim(wt);
  smwt <- (matrix(nrow=size[1], ncol=size[2]));
  smwt[size[1],] <- wt[size[1],]
  
  for (dt in (size[1] - 1):1 ){
    sm_table1 = sapply(pfOut1[dt,],function(x) dnorm(pfOut1[dt+1,],phi1 *  x,sqrt(1 - phi1^2)))
    bunsi = t(smwt[dt+1,] * sm_table1);
    bunbo = wt[dt,] %*% t(sm_table1);
    
    smwt[dt,] <- c(t(wt[dt,] * (bunsi %*% t(1/bunbo))))
  }
  return(smwt)
}


pairwise_weight <- function(phi1, filter_X, filter_weight, sm_weight){
  size = dim(filter_weight);
  pw_weight <- array(dim = c(size[2], size[2], size[1]));
  for (dt in 2:(size[1] - 1)){
    pw_table = sapply(filter_X[dt-1,],function(x) dnorm(filter_X[dt,],phi1 * x ,sqrt(1 - phi1^2)))
    bunbo = filter_weight[dt - 1,] %*% pw_table
    pw_weight[,,dt] <- t(t(pw_table * filter_weight[dt-1,]) * sm_weight[dt,] / c(bunbo))
  }
  return(pw_weight)
}


Q_calc <- function(par, pfOut1, rho1, pw_weight, sm_weight, y, v){
  phi1 <- sig(par[1])
  gam  <- exp(par[2])
  mu_g <- par[3]
  mu_f <- par[4]
  rho_f <- sig(par[5])
  V <- exp(par[6]);
  
  Q_state = 0;
  Q_obeserve = 0;
  first_state = 0;
  
  size <- dim(pfOut1)
  
  for (dt in 2:(size[1] - 1)){
    Q_state = Q_state + sum(sum(
      pw_weight[,,dt] * log(
        t(sapply(pfOut1[dt,],function(x) dnorm(x, phi1 * pfOut1[dt-1,], sqrt(1 - phi1^2))))
      )
    ));
    
    tmp1 = sapply(rho1[dt-1,],function(x) d_conditional_WJ(y[dt], y[dt-1], mu_g, x, mu_f, rho_f, 1));
    tmp2 = dgamma(v[dt]/(gam*exp(pfOut1[dt,]/2)) , shape = V, rate = V)/(gam*exp(pfOut1[dt,]/2));
    Q_obeserve = Q_obeserve + sm_weight[dt-1,] %*% log(tmp1) + sm_weight[dt,] %*% log(tmp2);
    
    if(is.nan(Q_state)|is.nan(Q_obeserve)){
      print(dt)
      browser()
    }
  }
  return(-Q_state - Q_obeserve);
}

Q_calc_para <- function(par, pfOut1, rho1, pw_weight, sm_weight, y, v){
  phi1 <- sig(par[1])
  gam  <- exp(par[2])
  mu_g <- par[3]
  mu_f <- par[4]
  rho_f <- sig(par[5])
  V <- exp(par[6]);
  mu_rho <- par[7];
  sig_rho <- exp(par[8]);
  
  Q_state = 0;
  Q_obeserve = 0;
  first_state = 0;
  
  size <- dim(pfOut1)
  Q_state = rep(0,size[1])
  
  clusterExport(scl,c('pw_weight','pfOut1','smwt','phi1','y','v','gam','V','mu_g','mu_f','rho_f','d_conditional_WJ',
                      'psswrappedcauchy','dsswrpcauchy','sig_rho','mu_rho'))
  
  tmp = parLapply(cl = scl,x = 2:(size[1] - 1),fun= function(dt){
    rho1_t_1<- 0.95 * ( tanh( sig_rho * pfOut1[dt-1,] + mu_rho)+1) / 2
    
    Q_state = sum(
      pw_weight[,,dt] * log(
        t(sapply(pfOut1[dt,],function(x) dnorm(x, phi1 * pfOut1[dt-1,], sqrt(1 - phi1^2))))
      )
    )
    
    tmp1 = sapply(rho1_t_1,function(x) d_conditional_WJ(y[dt], y[dt-1], mu_g, x, mu_f, rho_f, 1));
    tmp2 = dgamma(v[dt]/(gam*exp(pfOut1[dt,]/2)) , shape = V, rate = V)/(gam*exp(pfOut1[dt,]/2));
    Q_obeserve = Q_obeserve + sm_weight[dt-1,] %*% log(tmp1) + sm_weight[dt,] %*% log(tmp2);
    Q_state + Q_obeserve
  })
  tmp2 <- sapply(tmp,function(x)x)
  return(-sum(tmp2));
}
A <- matrix(c(1,2,3,4),ncol=2)
B <- c(1,3)



Q_calc_Estep <- function(par, pfOut1, rho1, pw_weight, sm_weight, y, v){
  phi1 <- par[1]
  gam  <- par[2]
  mu_g <- par[3]
  mu_f <- par[4]
  rho_f <- par[5]
  V <- par[6];
  mu_rho <- par[7];
  sig_rho <- par[8];
  
  Q_state = 0;
  Q_obeserve = 0;
  first_state = 0;
  
  size <- dim(pfOut1)
  
  for (dt in 2:(size[1] - 1)){
    rho1_t_1<- 0.95 * ( tanh( sig_rho * pfOut1[dt-1,] + mu_rho)+1) / 2
    Q_state = Q_state + sum(sum(
      pw_weight[,,dt] * log(
        t(sapply(pfOut1[dt,],function(x) dnorm(x, phi1 * pfOut1[dt-1,], sqrt(1 - phi1^2))))
      )
    ));
    
    tmp1 = sapply(rho1_t_1,function(x) d_conditional_WJ(y[dt], y[dt-1], mu_g, x, mu_f, rho_f, 1));
    tmp2 = dgamma(v[dt]/(gam*exp(pfOut1[dt,]/2)) , shape = V, rate = V)/(gam*exp(pfOut1[dt,]/2));
    Q_obeserve = Q_obeserve + sm_weight[dt-1,] %*% log(tmp1) + sm_weight[dt,] %*% log(tmp2);
    
    if(is.nan(Q_state)|is.nan(Q_obeserve)){
      print(dt)
      browser()
    }
  }
  return(-Q_state - Q_obeserve);
}