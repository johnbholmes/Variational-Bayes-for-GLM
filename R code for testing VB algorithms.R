#This R script contains code used to performing simulations in 
#An unified approach for Variational Bayesian inference in Generalised Linear Models by Holmes and Schofield.


####################################################################
## Starting values used
# 
# For all methods, starting values for $\hat{\boldsymbol{\beta}}$ and $\hat{\boldsymbol{\Sigma}}$ are chosen by assuming an intercept only model with estimate for $\hat{\boldsymbol{\beta}}_0$ chosen by applying 
# 
#
#In Poisson regression, the log transform to the sample mean of the data.
#In logistic regression, the logit  transform to the sample mean of the data.
#In negative binomial regression, we first apply the log  transform to the sample mean of the data. Second we calculate the probability of success by applying the inverse logit transform to $\log(E(Y)) - \log(k)$. 
#
########################### 
#Function inputs/outputs
# 
# For all methods the following inputs will be needed for all priors.
# 
# $\bf y$, the vector of responses of length $N$.
#$\bf X$, the predictor matrix (which includes intercept) of size $N \times p$.
#For logistic and negative binomial regression only, \texttt{Nsize}, the vector of binomial trial sizes/negative binomial failures of length $N$
#
# ###############
# Remaining function inputs are prior specific. These are 
# 
# Normal priors: $\boldsymbol{\beta}_p$, the vector of length $p$ of prior means and $\boldsymbol{\Sigma}_p^{-1}$, the prior inverse variance-covariance matrix. We have deliberately written the input as $\boldsymbol{\Sigma}_p^{-1}$ rather than $\boldsymbol{\Sigma}_p$, because a flat prior (likelihood alone in our setting) can then dealt with by setting $\boldsymbol{\Sigma}_p^{-1} = {\bf 0}_{p \times p}$.
# $t$ priors: $\boldsymbol{\mu}$, the $p$-length vector of prior degrees of freedom and ${\bf s}$, the $p$-length vector of prior scale parameters.
# Laplace priors: $\boldsymbol{\lambda}$, the $p$-length vector of laplace scale parameters. 
# 
# ####################################
# Outputs common to all methods are:
#
# $\hat{\boldsymbol{\beta}}$, the mean vector of the normal approximate posterior for $Q(\boldsymbol{\beta})$,
#  $\hat{\boldsymbol{\Sigma}}$, the variance-covariance matrix of the normal approximate posterior for $Q(\boldsymbol{\beta})$ 
# max_iter, the number of iterations required to reach convergence.
#
# ##################
# For methods where augmented variables are used, (MFVB, FF-MFVB), the additional outputs corresponding to the augmented variables are provided.
# $\hat{\boldsymbol{\omega}}$, the mean vector of the P{\'o}lya-Gamma augmented variables needed to implement MFVB versions of logistic and negative binomial regression.
# $\hat{\boldsymbol{\tau}}$, the $p$-length vector of expected precisions that appear $t$ and Laplace priors are represented as normal-gamma and normal-exponential mixtures respectively.
# 


#############################################################
#Background functions.

#To evaluate logit-normal moments, a C++ function has been developed: logitnormalmomentsderivatives1.cpp
#This can be found on the Github site: https://github.com/johnbholmes/Variational-Bayes-for-GLM, along with this R script.

#To evaluate the numerical integral needed to work directly with t-priors. The following R functions are required.

#integrate functions needed for VB t prior.
Ebt1<-function(input1,input2,input3){
  result<-integrate(f= function(x){x/(input3+x^2)*dnorm(x,mean=input1,sd=input2)},lower=input1-10*input2,upper=input1+10*input2)  
  return(result$val)  	}

Ebt2<-function(input1,input2,input3){
  result<-integrate(f= function(x){(input3-x^2)/(input3+x^2)^2*dnorm(x,mean=input1,sd=input2)},
                    lower=input1-10*input2,upper=input1+10*input2)  
  return(result$val)  	}




##############################Maximum a posteriori

####################
#Poisson regression

#Normal and flat prior		
Poisson.NMAP <-function(y,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-myp
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EY) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EY))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}
  }
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors
Poisson.tMAP <-function(y,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-myp
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}
  }
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Poisson.LMAP <-function(y,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-myp
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}
  }
  
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#############################
#Logistic regression

#Normal and flat prior
Logistic.NMAP <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-Nsize*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

# t-priors
Logistic.tMAP <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-Nsize*myp*(1-myp)
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}
  }
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Logistic.LMAP <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-Nsize*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}	

###############
# Negative binomial regression

#Normal and flat priors
Negbin.NMAP <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1)
  W0<-(y+Nsize)*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors
Negbin.tMAP <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1)
  W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Negbin.LMAP <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1)
  W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}	


##############################Fixed form Variational Bayes

#Note: This is fixed form variational Bayes, where we do not find MAP estimates first. This is slower than performing FFVB after finding MAP.

###############
#Poisson regression


#Normal and flat priors
Poisson.NVB <-function(y,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;W0<-myp
  priorFixed<-Sigmainvp%*%betap
  Sigma0     <- chol2inv(chol(crossprod(X*sqrt(W0))+Sigmainvp))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EY) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EY)) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors.
Poisson.tVB <-function(y,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;W0<-myp
  Sigma0     <- chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part #Note need functions Ebt1, Ebt2 loaded first. 
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    
    Fbeta<-crossprod(X,y-EY) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#Mimic t with logistic for faster computation
Poisson.tmimicVB <-function(y,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6
  myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters,observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;W0<-myp
  Sigma0     <- chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    
    Fbeta<-crossprod(X,y-EY) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Poisson.LVB <-function(y,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  myerr = rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0;W0<-myp
  Sigma0     <- chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part for logistic approximation of laplace
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior 
Poisson.trueLVB <-function(y,X, lambda){
  iter = 1000  #Maximum number of iterations
  myerr = rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0;W0<-myp
  Sigma0     <- chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part using true LASSO
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Square root of variance of beta needed in laplace penalty.
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-EY) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##########################
# Logistic regression


#Normal and flat priors
Logistic.NVB <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));	beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))+Sigmainvp))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,])) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#Independent t priors.
Logistic.tVB <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  
  #Iterations, VB part. Note, function Ebt1, Ebt2 need to be loaded first.
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#Mimic t with logistic for faster computation
Logistic.tmimicVB <-function(y,Nsize,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6
  myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Logistic.LVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter+1)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 	
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic for MAP only)
Logistic.trueLVB <-function(y,Nsize,X, lambda){
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))	
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Square root of variance of beta needed in laplace penalty.
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

###############################
#Negative binomial regression


Negbin.NVB <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))+Sigmainvp))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,])) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#Independent t priors.
Negbin.tVB <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1); W0<-(y+Nsize)*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part. Note functions Ebt1, Ebt2 need to be loaded first.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)		
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#Mimic t with logistic for faster computation
#Independent t priors.
Negbin.tmimicVB <-function(y,Nsize,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6
  myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}


##Laplace prior (LASSO mimic based on generalised logistic)
Negbin.LVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic for MAP only)
Negbin.trueLVB <-function(y,Nsize,X, lambda){	
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  Sigma0<-chol2inv(chol(crossprod(X*sqrt(W0))))
  
  #Iterations, VB part
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i+1]<-sum((beta0-betaold)^2)
    if( myerr[i+1] < (2/(p+3))*1e-6 | abs(myerr[i+1] - myerr[i]) < 1e-15) {break}	}
  if(myerr[i+1]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#########################################################################
#Fixed form Variational Bayes FFVB starting with MAP}
#Note: This is our preferred approach for fixed form variational Bayes in GLM. We start by finding MAP estimates first, to improve speed,
#particularly in logistic/negative binomial regression. In this section, we have included FFVB implementations of the type III logistic 
#which serve as approximations to the $t$ and Laplace priors. We recommend this for t priors, and the MAP part of the Laplace.


######################
#Poisson regression


#Normal and flat priors
Poisson.NMAPVB <-function(y,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0;	W0<-myp
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EY) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EY))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EY) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EY)) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors.
Poisson.tMAPVB <-function(y,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0;	W0<-myp
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    
    Fbeta<-crossprod(X,y-EY) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}







#Mimic t with logistic 
Poisson.tmimicMAPVB <-function(y,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6
  myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;	W0<-myp
  
  #Iterations, MAP part if logistic prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) + s*myscale - 2*s*myscale*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*s*myscale^2*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    Fbeta<-crossprod(X,y-EY) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Poisson.LMAPVB <-function(y,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  myerr = rep(0,iter)
  p<-dim(X)[2];n<-dim(X)[1]#no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);	beta0<-rep(0,p);beta0[1]<-b0;	W0<-myp
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}


##Laplace prior 
Poisson.trueLMAPVB <-function(y,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic for MAP. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  myerr = rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0;	W0<-myp
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )	
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part using true LASSO
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Square root of variance of beta needed in laplace penalty.
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-EY) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##################################
#Logistic regression

#Normal and flat priors
Logistic.NMAPVB <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations MAP parts.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,])) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors.
Logistic.tMAPVB <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0;	W0<-Nsize*myp*(1-myp)
  
  #Iterations, MAP part
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)		
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)		
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}





#Mimic t with logistic for faster computation
Logistic.tmimicMAPVB <-function(y,Nsize,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6
  myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));	beta0<-rep(0,p);beta0[1]<-b0;	W0<-Nsize*myp*(1-myp)
  
  #Iterations, MAP part if logistic prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + s*myscale - 2*s*myscale*exp(beta0*myscale)/(1+exp(beta0*myscale) )		
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*s*myscale^2*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 		
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Logistic.LMAPVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2];n<-dim(X)[1] #no parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n));b0 <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)#use for starting value
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )	
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 	
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior
Logistic.trueLMAPVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. This is needed for MAP. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  beta0<-rep(0,p);beta0[1]<-b0;	W0<-Nsize*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )		
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    #Square root of variance of beta needed in laplace penalty.
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)	}	

##############################
# Negative binomial regression

#Normal and flat priors
Negbin.NMAPVB <-function(y,Nsize,X, betap,Sigmainvp){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + priorFixed - Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,])) +Sigmainvp
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

#t priors.
Negbin.tMAPVB <-function(y,Nsize,X, nu,sj){
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Expectation of t prior.
    myEBt1 <-mapply(function(a,b,d){Ebt1(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)
    myEBt2 <-mapply(function(a,b,d){Ebt2(input1=a,input2=b,input3=d)},a=beta0,b=sqrt(diag(Sigma0)),d=nu*sj^2)		
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) -(nu+1)*myEBt1
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+(nu+1)*myEBt2
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}



#Mimic t with logistic 
Negbin.tmimicMAPVB <-function(y,Nsize,X, nu,sj){
  #Calculate approximate of t with logistic.
  r=2*(3*(nu+1))^0.5/nu;s=nu/6;myscale<-r/sj
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];	n<-dim(X)[1] #no parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations, MAP part if logistic prior.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + s*myscale - 2*s*myscale*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*s*myscale^2*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate t.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 	
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + s*myscale - 2*s*myscale*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*s*myscale^2*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior (LASSO mimic based on generalised logistic)
Negbin.LMAPVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05;myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2];	n<-dim(X)[1] #no parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )	
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  i0<-i
  #Iterations, VB part truly Laplace.
  for(i in (i0+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    bs <-beta0*myscale
    sbs<-sqrt(diag(Sigma0))*myscale
    EPbdiffPb <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=bs,b=sbs) 	
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) + lambda - 2*lambda*EPbdiffPb[1,]
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*myscale*EPbdiffPb[2,]
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##Laplace prior 
Negbin.trueLMAPVB <-function(y,Nsize,X, lambda){
  #Step 1: convert lambda to corresponding approximation
  #of logistic. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  myerr<-rep(0,iter)
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  myp  <-(1+exp(-b0+log(Nsize)))^(-1);	W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )	
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  i0<-i
  #Iterations, VB part
  for(i in (i0+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){
      EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Calculate moments of \beta in generalised logistic prior used to approximate laplace.
    sbs<-sqrt(diag(Sigma0))
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) +2*lambda*pnorm(-beta0/sbs) - lambda #penalised zero equation
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<- diag(Jbbeta)+ 2*lambda*dnorm(0,-beta0,sbs)
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0;Sigmaold<-Sigma0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    myerr[i]<-sum((beta0-betaold)^2)
    if( myerr[i] < (2/(p+3))*1e-6 | abs(myerr[i] - myerr[i-1]) < 1e-15) {break}	}
  if(myerr[i]> (2/(p+3))*1e-6){beta0<-0.5*(betaold+beta0);Sigma0<-0.5*(Sigma0+Sigmaold)}
  myresult<-list(beta=beta0,Sigma=Sigma0,max_iter=i)
  return(myresult)}

##############################################################################
# FF-MFVB starting with MAP

#In this section, we provide a mixed FFVB-MFVB implementation for GLM. This is only for $t$/Laplace priors and assumes the factorisation

#Poisson regression
#t priors
VBMAPt.MFFFVBPoisson<-function(X,y,nu,sj){
  epsilon = 1e-6
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-myp
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i
  tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
  
  #starting value for beta is MAP estimate.
  #Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}


#Laplace prior
VBMAPL.MFFFVBPoisson<-function(X,y,lambda){
  epsilon = 1e-6
  #Step 1: convert lambda to corresponding approximation
  #of logistic. This is needed for MAP only. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-mean(y) #use for starting value
  b0   <-log(myp)
  beta0<-rep(0,p);beta0[1]<-b0
  W0<-myp
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EY<-as.numeric(exp(Xb))
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EY)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i
  tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
  
  #starting value for beta is MAP estimate.
  #Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0
    xSx<-rowSums(X%*%Sigma0*X) #need sqrt of diagonal elements of X Sigma t(X).
    EY <-exp(as.numeric(Xb)+0.5*xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EY) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt(EY))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    #Update tau.
    tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}

#####################
# Logistic regression
#t prior
VBMAPt.MFFFVBlogistic<-function(X,y,Nsize,nu,sj){
  epsilon = 1e-6
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i;tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
  #starting value for beta is MAP estimate. Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}

#Laplace prior
VBMAPL.MFFFVBlogistic<-function(X,y,Nsize,lambda){
  epsilon = 1e-6
  #Step 1: convert lambda to corresponding approximation
  #of logistic. This is needed for MAP only. To stop things going too badly, assume alpha 0.05
  alpha=0.05
  myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no parameters/observations
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i;tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
  #starting value for beta is MAP estimate. Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){	EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-Nsize*EPdiffP[1,]) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt(Nsize*EPdiffP[2,]))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    #Update tau.
    tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}

##############
# Negative binomial regression
#t prior
VBMAPt.MFFFVBNegbin<-function(X,y,Nsize,nu,sj){
  epsilon = 1e-6
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameter/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i;	tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
  #starting value for beta is MAP estimate.Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}

#Laplace prior
VBMAPL.MFFFVBNegbin<-function(X,y,Nsize,lambda){
  epsilon = 1e-6
  #Step 1: convert lambda to corresponding approximation of logistic. Needed for MAP only. To stop things going too badly, assume alpha 0.05
  alpha=0.05;	myscale<-lambda/alpha #scale on generalised logistic
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2];n<-dim(X)[1] #no. parameters/observations
  myp<-mean(y) #use for starting value
  b0   <-log(myp);beta0<-rep(0,p);beta0[1]<-b0;myp  <-(1+exp(-b0+log(Nsize)))^(-1);W0<-(y+Nsize)*myp*(1-myp)
  
  #Iterations, MAP part if laplace (mimic) prior.
  for(i in 1:iter){
    Xb<-X%*%beta0-log(Nsize)
    EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
    Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i;tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
  #starting value for beta is MAP estimate. Iterations, VB part
  for(i in (imap+1):iter){
    Xb<-X%*%beta0-log(Nsize)
    xSx<-sqrt(rowSums(X%*%Sigma0*X)) #need sqrt of diagonal elements of X Sigma t(X).
    EPdiffP <-mapply(function(a,b){EPdiff01_logitnormalkpartcplus(mu=a,sigma=b)},a=Xb,b=xSx)
    #Setting up the zero equation and second derivative.
    Fbeta<-crossprod(X,y-(Nsize+y)*EPdiffP[1,]) - tau*beta0
    Jbbeta      <- crossprod(X*sqrt((Nsize+y)*EPdiffP[2,]))
    diag(Jbbeta)<-diag(Jbbeta)+tau
    #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    #Update tau.
    tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
    if( sum((beta0-betaold)^2) < (2/(p+3))*1e-6) {break} 	}
  
  param<-list(beta0,Sigma0,tau,i)
  names(param)<-c('beta','Sigma','tau','iter')
  return(param)}


##########################################################################################
####Mean field Variational Bayes (MFVB) starting with MAP}

#This is implementing  P{\'o}lya-Gamma data augmented Variational Bayes. This can only be used for logistic and negative
#binomial regression, as the method of \citet{Polson2013} is based on representing a type III logistic distribution kernel 
#as a normal-P{\'o}lya-Gamma mixture. The implied factorisation of the joint posterior in the MFVB approximation is:

#Q(\boldsymbol{\beta})\prod_{i=1}^NQ(\boldsymbol{\omega}_i)&& \quad \text{if prior is flat or Normal,} \\
#Q(\boldsymbol{\beta})\prod_{i=1}^NQ(\boldsymbol{\omega}_i)\prod_{j=1}^pQ(\tau_j)&& \quad \text{if prior is $t$ or Laplace.} 

#Details on the derivation are in the appendix of An unified approach for Variational Bayesian inference in Generalised Linear Models
#by Holmes and Schofield.

#Note: For logistic regression, this is a equivalent construction to the well known method of Jaakkola (1997).

###########################################################
# Logistic regression
#Normal and flat priors
VBMAPNorm.PGlogistic<-function(X,y,Nsize,betap,Sigmainvp){
  epsilon = 1e-6
  iter = 1000  #Maximum number of iterations
  p<-dim(X)[2] #no. parameters
  n<-dim(X)[1]
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp));beta0<-rep(0,p);beta0[1]<-b0;W0<-Nsize*myp*(1-myp)
  priorFixed<-Sigmainvp%*%betap #fixed part of prior
  
  #Iterations MAP parts.
  for(i in 1:iter){
    Xb<-X%*%beta0
    EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
    EdiffP<-EP*(Nsize-EP)/Nsize
    #Setting up the zero equation and second derivative.
    Sigmainvbeta <-Sigmainvp%*%beta0 
    Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
    Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
    betaold<-beta0
    Sigma0 <-chol2inv(chol(Jbbeta))
    beta0<-betaold + Sigma0%*%Fbeta
    if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
  
  imap<-i
  kappa <- y -0.5*Nsize #centred data for Jaakkola VB.
  
  #starting point for omega. 
  Xb <- X%*%beta0 
  VXb <- X%*%Sigma0%*%t(X)
  Exb <- sqrt(Xb^2 + diag(VXb))
  omega0 <- 0.5*Nsize*tanh(0.5*Exb)/Exb
  omega0 <- as.numeric(omega0)
  #Unchanged part for solving beta later. 
  XTY<-crossprod(X,kappa) + priorFixed
  XTX<-crossprod(X*sqrt(omega0)) + Sigmainvp #inverse variance-covariance.
  Vb0<-solve(XTX)
  
  #starting value for beta is MAP estimate.
  
  for(i in (imap+1):iter){
    b  <-Vb0%*%XTY     #update E(beta)
    b  <-as.numeric(b)
    #Update omega0
    Xb <- X%*%b 
    VXb <- X%*%Vb0%*%t(X)
    Exb <- sqrt(Xb^2 + diag(VXb))
    omega <- 0.5*Nsize*tanh(0.5*Exb)/Exb
    omega <-as.numeric(omega)
    Vb    <-chol2inv(chol(crossprod(X*sqrt(omega))+Sigmainvp ))
    
    diffb  <- sum((b-beta0)^2)
    if( diffb < (2/(p+3))*epsilon) break
    
    Vb0<-Vb;beta0<-b;omega0<-omega
    #Calculate relative change.	}
    
    param<-list(b,Vb,omega,i)
    names(param)<-c('beta','Sigma','omega','iter')
    return(param)}
  
  VBMAPt.PGlogistic<-function(X,y,Nsize,nu,sj){
    epsilon = 1e-6
    iter = 1000  #Maximum number of iterations
    p<-dim(X)[2] #no. parameters
    n<-dim(X)[1]
    myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
    b0   <-log(myp/(1-myp))
    beta0<-rep(0,p);beta0[1]<-b0
    W0<-Nsize*myp*(1-myp)
    
    #Iterations, MAP part
    for(i in 1:iter){
      Xb<-X%*%beta0
      EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
      EdiffP<-EP*(Nsize-EP)/Nsize
      #Setting up the zero equation and second derivative.
      Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
      
      Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
      diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
      betaold<-beta0
      Sigma0 <-chol2inv(chol(Jbbeta))
      beta0<-betaold + Sigma0%*%Fbeta
      if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
    
    imap<-i
    kappa <- y -0.5*Nsize #centred data for Jaakkola VB.
    
    #starting point for omega. 
    Xb <- X%*%beta0 
    VXb <- X%*%Sigma0%*%t(X)
    Exb <- sqrt(Xb^2 + diag(VXb))
    omega0 <- 0.5*Nsize*tanh(0.5*Exb)/Exb
    omega0 <- as.numeric(omega0)
    #Unchanged part for solving beta later. 
    XTY<-crossprod(X,kappa) 
    #Initial value for tau.
    tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
    #Initial value for Var(beta)
    XTX<-crossprod(X*sqrt(omega0)) #inverse variance-covariance.
    diag(XTX)<-diag(XTX)+tau
    Vb0<-solve(XTX)
    
    #starting value for beta is MAP estimate.
    
    for(i in (imap+1):iter){
      b  <-Vb0%*%XTY     #update E(beta)
      b  <-as.numeric(b)
      #Update omega0
      Xb <- X%*%b 
      VXb <- X%*%Vb0%*%t(X)
      Exb <- sqrt(Xb^2 + diag(VXb))
      omega <- 0.5*Nsize*tanh(0.5*Exb)/Exb
      omega <-as.numeric(omega)
      #Update Var(beta)
      VBinv<-crossprod(X*sqrt(omega))
      diag(VBinv)<-diag(VBinv) + tau
      Vb    <-chol2inv(chol(VBinv))
      #Update tau.
      tau  <- (nu+1)/(nu*sj^2+b^2+diag(Vb))
      
      diffb  <- sum((b-beta0)^2)
      if( diffb < (2/(p+3))*epsilon) break
      
      Vb0<-Vb;beta0<-b;omega0<-omega;tau0<-tau	}
    
    param<-list(b,Vb,omega,tau,i)
    names(param)<-c('beta','Sigma','omega','tau','iter')
    return(param)}
  
  #Laplace prior
  VBMAPL.PGlogistic<-function(X,y,Nsize,lambda){
    epsilon = 1e-6
    #Step 1: convert lambda to corresponding approximation
    #of logistic. To stop things going too badly, assume alpha 0.05
    alpha=0.05
    myscale<-lambda/alpha #scale on generalised logistic
    
    iter = 1000  #Maximum number of iterations
    p<-dim(X)[2] #no. parameters
    n<-dim(X)[1]
    myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
    b0   <-log(myp/(1-myp))
    beta0<-rep(0,p);beta0[1]<-b0
    W0<-Nsize*myp*(1-myp)
    
    #Iterations, MAP part if laplace (mimic) prior.
    for(i in 1:iter){
      Xb<-X%*%beta0
      EP<-Nsize*as.numeric((1+exp(-Xb))^(-1))
      EdiffP<-EP*(Nsize-EP)/Nsize
      #Setting up the zero equation and second derivative.
      Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
      
      Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
      diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
      betaold<-beta0
      Sigma0 <-chol2inv(chol(Jbbeta))
      beta0<-betaold + Sigma0%*%Fbeta
      if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
    
    imap<-i
    kappa <- y -0.5*Nsize #centred data for Jaakkola VB.
    
    #starting point for omega. 
    Xb <- X%*%beta0 
    VXb <- X%*%Sigma0%*%t(X)
    Exb <- sqrt(Xb^2 + diag(VXb))
    omega0 <- 0.5*Nsize*tanh(0.5*Exb)/Exb
    omega0 <- as.numeric(omega0)
    #Unchanged part for solving beta later. 
    XTY<-crossprod(X,kappa) 
    #Initial value for tau.
    tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
    #Initial value for Var(beta)
    XTX<-crossprod(X*sqrt(omega0)) #inverse variance-covariance.
    diag(XTX)<-diag(XTX)+tau
    Vb0<-solve(XTX)
    
    #starting value for beta is MAP estimate.
    for(i in (imap+1):iter){
      b  <-Vb0%*%XTY     #update E(beta)
      b  <-as.numeric(b)
      #Update omega0
      Xb <- X%*%b 
      VXb <- X%*%Vb0%*%t(X)
      Exb <- sqrt(Xb^2 + diag(VXb))
      omega <- 0.5*Nsize*tanh(0.5*Exb)/Exb
      omega <-as.numeric(omega)
      #Update Var(beta)
      VBinv<-crossprod(X*sqrt(omega))
      diag(VBinv)<-diag(VBinv) + tau
      Vb    <-chol2inv(chol(VBinv))
      #Update tau.
      tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
      
      diffb  <- sum((b-beta0)^2)
      if( diffb < (2/(p+3))*epsilon) break
      
      Vb0<-Vb;beta0<-b;omega0<-omega;tau0<-tau	}
    
    param<-list(b,Vb,omega,tau,i)
    names(param)<-c('beta','Sigma','omega','tau','iter')
    return(param)}


#######################  
#Negative binomial regression

  #Normal and flat prior		
  VBMAPNorm.PGNegbin<-function(X,y,Nsize,betap,Sigmainvp){
    epsilon = 1e-6
    iter = 1000  #Maximum number of iterations
    p<-dim(X)[2] #no. parameters
    n<-dim(X)[1]
    myp<-mean(y) #use for starting value
    b0   <-log(myp)
    beta0<-rep(0,p);beta0[1]<-b0
    myp  <-(1+exp(-b0+log(Nsize)))^(-1)
    W0<-(y+Nsize)*myp*(1-myp)
    priorFixed<-Sigmainvp%*%betap
    
    #Iterations
    for(i in 1:iter){
      Xb<-X%*%beta0-log(Nsize)
      EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
      EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
      #Setting up the zero equation and second derivative.
      Sigmainvbeta <-Sigmainvp%*%beta0 
      Fbeta<-crossprod(X,y-EP) + priorFixed -Sigmainvbeta
      Jbbeta      <- crossprod(X*sqrt(EdiffP))+Sigmainvp #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
      betaold<-beta0
      Sigma0 <-chol2inv(chol(Jbbeta))
      beta0<-betaold + Sigma0%*%Fbeta
      if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
    
    imap<-i
    kappa <- y -0.5*(Nsize+y) #centred data for Jaakkola VB.
    
    #starting point for omega. 
    Xb <- X%*%beta0 -log(Nsize)
    VXb <- X%*%Sigma0%*%t(X)
    Exb <- sqrt(Xb^2 + diag(VXb))
    omega0 <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
    omega0 <- as.numeric(omega0)
    #Unchanged part for solving beta later. 
    XTY<-crossprod(X,kappa) + priorFixed
    XTX<-crossprod(X*sqrt(omega0)) + Sigmainvp #inverse variance-covariance.
    Vb0<-solve(XTX)
    
    #starting value for beta is MAP estimate.
    XTY2<-XTY+crossprod(X,omega0*log(Nsize))
    for(i in (imap+1):iter){
      
      b  <-Vb0%*%XTY2     #update E(beta)
      b  <-as.numeric(b)
      #Update omega0
      Xb <- X%*%b -log(Nsize)
      VXb <- X%*%Vb0%*%t(X)
      Exb <- sqrt(Xb^2 + diag(VXb))
      omega <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
      omega <-as.numeric(omega)
      Vb    <-chol2inv(chol(crossprod(X*sqrt(omega))+Sigmainvp ))
      
      diffb  <- sum((b-beta0)^2)
      #update XTY
      XTY2 <- XTY +crossprod(X,omega*log(Nsize))
      if( diffb < (2/(p+3))*epsilon) break
      
      Vb0<-Vb;beta0<-b;omega0<-omega
      #Calculate relative change.	}
      
      param<-list(b,Vb,omega,i)
      names(param)<-c('beta','Sigma','omega','iter')
      return(param)}
    
    #t prior
    VBMAPt.PGNegbin<-function(X,y,Nsize,nu,sj){
      epsilon = 1e-6
      iter = 1000  #Maximum number of iterations
      p<-dim(X)[2] #no. parameters
      n<-dim(X)[1]
      myp<-mean(y) #use for starting value
      b0   <-log(myp)
      beta0<-rep(0,p);beta0[1]<-b0
      myp  <-(1+exp(-b0+log(Nsize)))^(-1)
      W0<-(y+Nsize)*myp*(1-myp)
      
      #Iterations
      for(i in 1:iter){
        Xb<-X%*%beta0-log(Nsize)
        EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
        EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
        #Setting up the zero equation and second derivative.
        Fbeta<-crossprod(X,y-EP) -(nu+1)*beta0/(nu*sj^2+beta0^2)
        
        Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
        diag(Jbbeta) <- diag(Jbbeta) + (nu+1)*(nu*sj^2-beta0^2)/(nu*sj^2+beta0^2)^2 #Include prior information.
        betaold<-beta0
        Sigma0 <-chol2inv(chol(Jbbeta))
        beta0<-betaold + Sigma0%*%Fbeta
        if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
      
      imap<-i
      kappa <- y -0.5*(Nsize+y) #centred data for Jaakkola VB.
      
      #starting point for omega. 
      Xb <- X%*%beta0 -log(Nsize)
      VXb <- X%*%Sigma0%*%t(X)
      Exb <- sqrt(Xb^2 + diag(VXb))
      omega0 <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
      omega0 <- as.numeric(omega0)
      #Unchanged part for solving beta later. 
      XTY<-crossprod(X,kappa) 
      #Initial value for tau.
      tau  <- (nu+1)/(nu*sj^2+beta0^2+diag(Sigma0))
      #Initial value for Var(beta)
      XTX<-crossprod(X*sqrt(omega0)) #inverse variance-covariance.
      diag(XTX)<-diag(XTX)+tau
      Vb0<-solve(XTX)
      
      XTY2<-XTY+crossprod(X,omega0*log(Nsize))
      #starting value for beta is MAP estimate.
      
      for(i in (imap+1):iter){
        b  <-Vb0%*%XTY2    #update E(beta)
        b  <-as.numeric(b)
        #Update omega0
        Xb <- X%*%b -log(Nsize)
        VXb <- X%*%Vb0%*%t(X)
        Exb <- sqrt(Xb^2 + diag(VXb))
        omega <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
        omega <-as.numeric(omega)
        #Update Var(beta)
        VBinv<-crossprod(X*sqrt(omega))
        diag(VBinv)<-diag(VBinv) + tau
        Vb    <-chol2inv(chol(VBinv))
        #Update tau.
        tau  <- (nu+1)/(nu*sj^2+b^2+diag(Vb))
        
        diffb  <- sum((b-beta0)^2)
        #update XTY
        XTY2 <- XTY +crossprod(X,omega*log(Nsize))
        if( diffb < (2/(p+3))*epsilon) break
        
        Vb0<-Vb;beta0<-b;omega0<-omega;tau0<-tau	}
      
      param<-list(b,Vb,omega,tau,i)
      names(param)<-c('beta','Sigma','omega','tau','iter')
      return(param)}
    
    #Laplace prior
    VBMAPL.PGNegbin<-function(X,y,Nsize,lambda){
      epsilon = 1e-6
      #Step 1: convert lambda to corresponding approximation
      #of logistic. To stop things going too badly, assume alpha 0.05
      alpha=0.05
      myscale<-lambda/alpha #scale on generalised logistic
      
      iter = 1000  #Maximum number of iterations
      p<-dim(X)[2] #no. parameters
      n<-dim(X)[1]
      myp<-mean(y) #use for starting value
      b0   <-log(myp)
      beta0<-rep(0,p);beta0[1]<-b0
      myp  <-(1+exp(-b0+log(Nsize)))^(-1)
      W0<-(y+Nsize)*myp*(1-myp)
      
      #Iterations, MAP part if laplace (mimic) prior.
      for(i in 1:iter){
        Xb<-X%*%beta0-log(Nsize)
        EP<-(Nsize+y)*as.numeric((1+exp(-Xb))^(-1))
        EdiffP<-EP*(Nsize+y-EP)/(Nsize+y)
        #Setting up the zero equation and second derivative.
        Fbeta<-crossprod(X,y-EP) + lambda - 2*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )
        
        Jbbeta       <- crossprod(X*sqrt(EdiffP)) #Finding t(X)diag(W)X= t(X Hadamard mult W)X.
        diag(Jbbeta) <- diag(Jbbeta) + 2*myscale*lambda*exp(beta0*myscale)/(1+exp(beta0*myscale) )^2 #Include prior information.
        betaold<-beta0
        Sigma0 <-chol2inv(chol(Jbbeta))
        beta0<-betaold + Sigma0%*%Fbeta
        if( sum((beta0-betaold)^2) < (2/(length(beta0)+3))*1e-6) {break}	}
      
      imap<-i
      kappa <- y -0.5*(Nsize+y) #centred data for Jaakkola VB.
      
      #starting point for omega. 
      Xb <- X%*%beta0 -log(Nsize)
      VXb <- X%*%Sigma0%*%t(X)
      Exb <- sqrt(Xb^2 + diag(VXb))
      omega0 <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
      omega0 <- as.numeric(omega0)
      #Unchanged part for solving beta later. 
      XTY<-crossprod(X,kappa) 
      #Initial value for tau.
      tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
      #Initial value for Var(beta)
      XTX<-crossprod(X*sqrt(omega0)) #inverse variance-covariance.
      diag(XTX)<-diag(XTX)+tau
      Vb0<-solve(XTX)
      
      XTY2<-XTY+crossprod(X,omega0*log(Nsize))
      #starting value for beta is MAP estimate.
      
      for(i in (imap+1):iter){
        b  <-Vb0%*%XTY2     #update E(beta)
        b  <-as.numeric(b)
        #Update omega0
        Xb <- X%*%b -log(Nsize)
        VXb <- X%*%Vb0%*%t(X)
        Exb <- sqrt(Xb^2 + diag(VXb))
        omega <- 0.5*(Nsize+y)*tanh(0.5*Exb)/Exb
        omega <-as.numeric(omega)
        #Update Var(beta)
        VBinv<-crossprod(X*sqrt(omega))
        diag(VBinv)<-diag(VBinv) + tau
        Vb    <-chol2inv(chol(VBinv))
        #Update tau.
        tau  <- lambda/sqrt(beta0^2+diag(Sigma0))
        
        diffb  <- sum((b-beta0)^2)
        #update XTY
        XTY2 <- XTY +crossprod(X,omega*log(Nsize))
        if( diffb < (2/(p+3))*epsilon) break
        
        Vb0<-Vb;beta0<-b;omega0<-omega;tau0<-tau	}
      
      param<-list(b,Vb,omega,tau,i)
      names(param)<-c('beta','Sigma','omega','tau','iter')
      return(param)}

