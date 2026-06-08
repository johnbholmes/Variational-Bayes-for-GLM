#######################################################################################################
# Polya-Gamma data augmented logistic and negative binomial regression

#Inputs for PG logit/negbin: 
#y number of successes
#Nsize number of trial (logistic) or number of failures (negative binomial)
#X predictor matrix
#iter, no. of iterations
#burnin, initial iterations to throw out.
#family: Two options 'logistic' or 'negative binomial'
#Now run PG models for the various models
#This will cover both logistic and negative binomial models

#for flat  priors.
PGflat.logit<-function(y,Nsize,X,iter,burnin,family){
  Negbin.ind<-family != 'logistic' #indicator for negative binomial or logistic
  #1 if neg-bin, 0 if logistic
  #Offset if negative binomial
  logN   <-Negbin.ind*log(Nsize)
  #'Totals' if negative binomial 
  Nsize  <-Nsize + Negbin.ind*y #modify Nsize if Negative binomial

  library(BayesLogit)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p+n)
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n))*(1-Negbin.ind) +mean(y)*Negbin.ind#use for starting value
  b0   <-log(myp/(1-myp)*(1-Negbin.ind) + myp*Negbin.ind)
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    #If negative binomial need to add omega*logN, if logistic 
    #logN is set to zero.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega)+omega*logN) 
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X)
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb-logN) 
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega)  
  }
  return(betavarcomp.store)  
}

PGNorm.logit<-function(y,Nsize,X,iter,burnin,betap,Sigmainvp,family){
  Negbin.ind<-family != 'logistic' #indicator for negative binomial or logistic
  #1 if neg-bin, 0 if logistic
  #Offset if negative binomial
  logN   <-Negbin.ind*log(Nsize)
  #'Totals' if negative binomial 
  Nsize  <-Nsize + Negbin.ind*y #modify Nsize if Negative binomial
  
  library(BayesLogit)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  priorFixed<-Sigmainvp%*%betap
  pbeta <- pbeta+priorFixed
  #Cholesky of prior inverse for simulating noise.
  Sigmainvchol <- chol(Sigmainvp)
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p+n)
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n))*(1-Negbin.ind) +mean(y)*Negbin.ind#use for starting value
  b0   <-log(myp/(1-myp)*(1-Negbin.ind) + myp*Negbin.ind)
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega)+omega*logN)  + t(Sigmainvchol)%*%rnorm(p)
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X) + Sigmainvp
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb-logN) 
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega)  
  }
  return(betavarcomp.store)  
}

PGt.logit<-function(y,Nsize,X,iter,burnin,nu,sj,family){
  Negbin.ind<-family != 'logistic' #indicator for negative binomial or logistic
  #1 if neg-bin, 0 if logistic
  #Offset if negative binomial
  logN   <-Negbin.ind*log(Nsize)
  #'Totals' if negative binomial 
  Nsize  <-Nsize + Negbin.ind*y #modify Nsize if Negative binomial
  
  library(BayesLogit)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  #draw tau from prior.
  tau  <-rgamma(p,0.5*(nu+1),0.5*(nu*sj^2))
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p*2+n) #as we have augmented variables 
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n))*(1-Negbin.ind) +mean(y)*Negbin.ind#use for starting value
  b0   <-log(myp/(1-myp)*(1-Negbin.ind) + myp*Negbin.ind)
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega)+omega*logN)  + rnorm(p)*sqrt(tau)
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X) 
    diag(Sigmainv)<-diag(Sigmainv)+tau
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb-logN) 
    tau  <-rgamma(p,0.5*(nu+1),0.5*(nu*sj^2+beta^2))
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega,tau)  
  }
  return(betavarcomp.store)  
}

#Bayesian Logistic LASSO, fixed lambda
PGL.logit<-function(y,Nsize,X,iter,burnin,lambda,family){
  Negbin.ind<-family != 'logistic' #indicator for negative binomial or logistic
  #1 if neg-bin, 0 if logistic
  #Offset if negative binomial
  logN   <-Negbin.ind*log(Nsize)
  #'Totals' if negative binomial 
  Nsize  <-Nsize + Negbin.ind*y #modify Nsize if Negative binomial
  
  library(BayesLogit)
  library(LaplacesDemon)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  #draw tau from prior.
  tau  <-rinvgaussian(p,lambda/abs(rnorm(p)),lambda^2) 
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p*2+n) #as we have augmented variables 
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n))*(1-Negbin.ind) +mean(y)*Negbin.ind#use for starting value
  b0   <-log(myp/(1-myp)*(1-Negbin.ind) + myp*Negbin.ind)
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega)+omega*logN)  + rnorm(p)*sqrt(tau)
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X) 
    diag(Sigmainv)<-diag(Sigmainv)+tau
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb-logN) 
    tau  <-rinvgaussian(p,lambda/abs(beta),lambda^2)
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega,1/tau)  
  }
  return(betavarcomp.store)  
}

#Gibbs sampler for the lasso mimic using generalised logistic prior.
PGLmic.logit<-function(y,Nsize,X,iter,burnin,lambda,family){
  Negbin.ind<-family != 'logistic' #indicator for negative binomial or logistic
  #1 if neg-bin, 0 if logistic
  #Offset if negative binomial
  logN   <-Negbin.ind*log(Nsize)
  #'Totals' if negative binomial 
  Nsize  <-Nsize + Negbin.ind*y #modify Nsize if Negative binomial
  
  #Step 1: convert lambda to corresponding approximation
  #of logistic.
  mys<- 0.05 
  myscale<-lambda/mys #scale on generalised logistic
  
  library(BayesLogit)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  #draw tau from prior. PG(2mys,0)
  tau  <-rpg(p,2*mys,0)
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p*2+n) #as we have augmented variables 
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n))*(1-Negbin.ind) +mean(y)*Negbin.ind#use for starting value
  b0   <-log(myp/(1-myp))*(1-Negbin.ind) + log(myp)*Negbin.ind
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega)+omega*logN)  + rnorm(p)*sqrt(tau)*myscale
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X) 
    diag(Sigmainv)<-diag(Sigmainv)+tau*myscale^2
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb-logN) 
    tau  <-rpg(p,2*mys,(myscale*beta))
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega,tau)  
  }
  return(betavarcomp.store)  
}
