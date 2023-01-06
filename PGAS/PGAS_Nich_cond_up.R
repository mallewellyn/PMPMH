########### load packages 
library(invgamma)

######### simulate data 
set.seed(66)
n=300
N=numeric(n)
R=numeric(n)
S=numeric(n)
y=numeric(n)
tau=5
P=50
N0=1000
delta=0.7
phi=1
sigma2E=10
sigma2Epsilon=1
betaE=1/sigma2E
betaEpsilon=1/sigma2Epsilon
e=rgamma(n, betaE, rate=betaE)
epsilon=rgamma(n, betaEpsilon, rate=betaEpsilon)


S[1]=rbinom(1, N0, exp(-delta*epsilon[1]))
N[1]=S[1]
for(t in 2:(tau)){
  S[t]=rbinom(1, N[t-1], exp(-delta*epsilon[t]))
  N[t]=S[t]
}
S[tau+1]=rbinom(1, N[tau], exp(-delta*epsilon[tau+1]))
R[tau+1]=rpois(1, P*N0*e[tau+1]*exp(-1))
N[tau+1]=S[tau+1]+R[tau+1]
for(i in (tau+2):n){
  R[i]=rpois(1, P*N[i-tau-1]*e[i]*exp(-N[i-tau-1]/N0))
  S[i]=rbinom(1, N[i-1], exp(-delta*epsilon[i]))
  N[i]=R[i]+S[i]
}
for(i in 1:n){
  y[i]=rpois(1, phi*N[i])
}

set.seed(NULL)


########### PGAS function
uniformSpacings<-function(N){
  S=numeric(N+1)
  E=numeric(N+1)
  for(n in 1:(N+1)){
    E[n]=rexp(1, 1)
    if(n==1){
      S[n]=E[n]
    } else {
      S[n]=S[n-1]+E[n]
    }
  }
  U=numeric(N)
  for(n in 1:N){
    U[n]=S[n]/S[N+1]
  }
  return(U)
}

icdf<-function(weights, U){
  N=length(U)
  A=numeric(N)
  s=weights[1]
  m=1
  
  for(n in 1:N){
    while(s<U[n]){
      m=m+1
      s=s+weights[m]
    }
    A[n]=m
  }
  return(A)
  
}

multinomialResample<-function(weights){
  N=length(weights)-1
  U<-uniformSpacings(N)
  A=icdf(weights, U)
  return(A)
}

PGASR<-function(N, rESSMin, y, theta, e, epsilon, RInit, SInit){
  
  particles<-matrix(rep(0, N*length(y)), ncol=length(y))
  weights<-matrix(rep(0, N*length(y)), ncol=length(y))
  A<-matrix(rep(0, N*length(y)), ncol=length(y))
  normWeights<-matrix(rep(0, N*length(y)), ncol=length(y))
  aWeights<-matrix(rep(0, N*length(y)), ncol=length(y))
  aNormWeights<-matrix(rep(0, N*length(y)), ncol=length(y))

  priorRandInit<-function(t, theta){
    rpois(1, theta[2]*N0*exp(-1)*e[6])
  }
  priorRandOtherwise<-function(t, theta, partMinusTau){
    NtMinusTau=SInit[t-tau-1] + partMinusTau
    rpois(1, theta[2]*NtMinusTau*exp(-NtMinusTau/N0)*e[t])
  }
  logLikelihood<-function(t, theta, partst){
    dens=numeric(length(partst))
    N=partst+SInit[t]
    if(t!=length(y)){
      for(i in 1:length(partst)){
        dens[i]=dbinom(SInit[t+1], N[i], exp(-theta[1]*epsilon[t+1]), log=TRUE) +
          dpois(y[t], theta[5]*(N[i]), log=TRUE)
      }
    } else {
      for(i in 1:length(partst)){
        dens[i]=dpois(y[t], theta[5]*N[i], log=TRUE)
      }
    }
    dens
    
  }

  
  weights[,1:tau]=1
  normWeights[,1:tau]=1/N
  #time tau+1
  particles[2:N,tau+1]<-sapply(rep(1, N-1), priorRandInit, theta=theta)
  particles[1,tau+1]=RInit[tau+1]
  
  logWeights<-logLikelihood(tau+1, theta, particles[,tau+1])
  if(length(which(logWeights==-Inf))==N){
    logWeights=rep(log(1), N)
  }
  const=max(logWeights)
  weights[,tau+1]=exp(logWeights-const)
  normWeights[,tau+1]=weights[,tau+1]/sum(weights[,tau+1])
  
  for(t in (tau+2):length(y)){
    if(((1/sum(normWeights[,t-1]^2))/N)<rESSMin){
      A[2:N,t]<-multinomialResample(normWeights[,t-1])
      weights[,t-1]=1
    } else {
      A[2:N,t]=2:N
    }
    
    particles[2:N,t]=sapply(particles[,t-tau-1][A[-1,t]], priorRandOtherwise, theta=theta, t=t)
    
    particles[1,t]=RInit[t]
    
    ##ancestor sampling 
    inds=seq(t, t+5, by=1)
    inds=inds[which(inds<=length(y))]
    
    RFunc<-function(l, i){
      dpois(RInit[l], theta[2]*e[l]*(particles[i,l-tau-1]+SInit[l-tau-1])*exp(-(particles[i,l-tau-1]+SInit[l-tau-1])/N0), log=T)
    }
    anFunc<-function(i){
      sum(sapply(inds, RFunc, i=i)) + dbinom(SInit[t], particles[i,t-1]+SInit[t-1], exp(-theta[1]*epsilon[t]), log=T)
    }
    aLogWeights=logWeights+sapply(1:N, anFunc)
    
    if(length(which(aLogWeights==-Inf))==N){
      aLogWeights=rep(log(1/N), N)
    }
    const=max(aLogWeights)
    aWeights[,t]<-exp(aLogWeights-const)
    aNormWeights[,t]<-aWeights[,t]/sum(aWeights[,t])
    A[1,t]=sample(multinomialResample(aNormWeights[,t]), size=1)
    
    ####
    particles[,1:t]=cbind(particles[A[,t],1:(t-1)], particles[,t])
    
    logWeights<-log(weights[A[,t],t-1])+ logLikelihood(t, theta, particles[,t])
    if(length(which(logWeights==-Inf))==N){
      logWeights=rep(log(1), N)
    }
    
    const=max(logWeights)
    weights[,t]<-exp(logWeights-const)
    normWeights[,t]<-weights[,t]/sum(weights[,t])
    
  }
  
  samp<-numeric(length(y))
  xx<-sample(multinomialResample(normWeights[,length(y)]), size=1)
  samp=particles[xx,]
  
  samp
}

PGASS<-function(N, rESSMin, y, theta, e, epsilon, SInit, RInit){
  
  particles<-matrix(rep(0, N*length(y)), ncol=length(y))
  weights<-matrix(rep(0, N*length(y)), ncol=length(y))
  A<-matrix(rep(0, N*length(y)), ncol=length(y))
  normWeights<-matrix(rep(0, N*length(y)), ncol=length(y))
  aWeights<-matrix(rep(0, N*length(y)), ncol=length(y))
  aNormWeights<-matrix(rep(0, N*length(y)), ncol=length(y))

  priorRandInit<-function(t, theta){
    rbinom(1, N0, exp(-theta[1]*epsilon[1]))
  }
  priorRandOtherwise<-function(t, theta, partMinus){
    NtMinus=RInit[t-1] + partMinus
    rbinom(1, NtMinus, exp(-theta[1]*epsilon[t]))
  }
  
  logLikelihood<-function(t, theta, partst){
    dens=numeric(length(partst))
    Nt=partst+RInit[t]
    
    if(t > length(y)-(tau+1)){
      for(i in 1:length(partst)){
        dens[i]=dpois(y[t], theta[5]*Nt[i], log=TRUE)
      }
    } else { 
      for(i in 1:length(partst)){
        dens[i]=dpois(y[t], theta[5]*Nt[i], log=TRUE) + 
          dpois(RInit[t+tau+1], theta[2]*Nt[i]*e[t+tau+1]*exp(-Nt[i]/N0), log=TRUE)
        
      }
    } 
    
    
    dens
    
  }
  
  #from time 1 this time
  particles[2:N,1]<-sapply(rep(1, N-1), priorRandInit, theta=theta)
  particles[1,1]=SInit[1]
  
  logWeights<-logLikelihood(1, theta, particles[,1])
  if(length(which(logWeights==-Inf))==N){
    logWeights=rep(log(1), N)
  }
  const=max(logWeights)
  weights[,1]=exp(logWeights-const)
  normWeights[,1]=weights[,1]/sum(weights[,1])
  
  for(t in 2:length(y)){
    if(((1/sum(normWeights[,t-1]^2))/N)<rESSMin){
      A[2:N,t]<-multinomialResample(normWeights[,t-1])
      weights[,t-1]=1
    } else {
      A[2:N,t]=2:N
    }
    
    particles[2:N,t]=sapply(particles[,t-1][A[-1,t]], priorRandOtherwise, theta=theta, t=t)
    
    particles[1,t]=SInit[t]
    
    ##ancestor sampling
    anFunc<-function(i){
      dbinom(SInit[t], particles[i,t-1] + RInit[t-1], exp(-theta[1]*epsilon[t]), log=T)
    }
    aLogWeights=logWeights+sapply(1:N, anFunc)
    
    if(length(which(aLogWeights==-Inf))==N){
      xx<-which(RInit[t-1] + particles[,t-1]>=SInit[t])
      aLogWeights=rep(log(1), N)
      aLogWeights[xx]=rep(log(1/length(xx)), length(xx))
    }
    const=max(aLogWeights)
    aWeights[,t]<-exp(aLogWeights-const)
    aNormWeights[,t]<-aWeights[,t]/sum(aWeights[,t])
    A[1,t]=sample(multinomialResample(aNormWeights[,t]), size=1)
    
    ##
    particles[,1:t]=cbind(particles[A[,t],1:(t-1)], particles[,t])
    
    logWeights<-log(weights[A[,t],t-1])+logLikelihood(t, theta, particles[,t])
    if(length(which(logWeights==-Inf))==N){
      logWeights=rep(log(1), N)
    }
    
    const=max(logWeights)
    weights[,t]<-exp(logWeights-const)
    normWeights[,t]<-weights[,t]/sum(weights[,t])
    
  }
  
  samp<-numeric(length(y))
  xx<-sample(multinomialResample(normWeights[,length(y)]), size=1)
  samp=particles[xx,]
  
  samp
}


########### Model parameter sampling 
deltaA=0.007
deltaB=0.01
PA=50
PB=1
betaEpsilonA=100
betaEpsilonB=100
betaEA=10
betaEB=1
phiA=0.01
phiB=0.01
propLim=c(0.03, NA, 0.5, 0.05, NA)

thetaUpdate<-function(SStates, RStates, theta){
  
  NCurr=SStates+RStates
  thetaCurr<-theta
  thetaProp<-theta
  pacc<-numeric(length(theta))
  
  logPriorDensityBetaEpsilon<-function(val, betaEpsilonA, betaEpsilonB){
    dinvgamma(val, shape=betaEpsilonA, rate=betaEpsilonB, log=TRUE)
  }
  
  logPriorDensityBetaE<-function(val, betaEA, betaEB){
    dinvgamma(val, shape=betaEA, rate=betaEB, log=TRUE)
  }
  
  ##delta
  thetaProp[1]=runif(1, thetaCurr[1]-propLim[1], thetaCurr[1] + propLim[1])
  
  logLik<-function(t, delta){
    if(t==1){
      dbinom(SStates[t], N0, exp(-delta*epsilon[t]), log=TRUE)
    } else {
      dbinom(SStates[t], NCurr[t-1], exp(-delta*epsilon[t]), log=TRUE)
    }
  }
  
  if(thetaProp[1]>0){
    pacc[1]=min(1, exp(dgamma(thetaProp[1], deltaA, rate=deltaB, log=TRUE) +
                         sum(sapply(1:length(y), logLik, delta=thetaProp[1])) - (
                           dgamma(thetaCurr[1], deltaA, rate=deltaB, log=TRUE) +
                             sum(sapply(1:length(y), logLik, delta=thetaCurr[1])) )))
  }
  if(is.na(pacc[1])){
    pacc[1]=0
  }
  if(runif(1)<pacc[1]){
    thetaCurr[1]=thetaProp[1]
  } else {
    thetaProp[1]=thetaCurr[1]
  }
  
  ## P
  thetaProp[2]=rgamma(1, PA + sum(RStates[(tau+1):length(y)]),
                      rate=PA + N0*exp(-1)*e[tau+1] + sum(NCurr[1:(length(y)-tau-1)]*exp(-NCurr[1:(length(y)-tau-1)]/N0)*e[(tau+2):length(y)]))
  thetaCurr[2]=thetaProp[2]
  pacc[2]=1
  
  #beta_epsilon
  thetaProp[3]=runif(1, thetaCurr[3]-propLim[3], thetaCurr[3]+propLim[3])
  if(thetaProp[3]>0){
    pacc[3]=min(1, exp(dinvgamma(thetaProp[3], betaEpsilonA, rate=betaEpsilonB, log=TRUE) + sum(dgamma(epsilon, shape=thetaProp[3], rate=thetaProp[3], log=TRUE)) -
                         (dinvgamma(thetaCurr[3], betaEpsilonA, rate=betaEpsilonB, log=TRUE) + sum(dgamma(epsilon, shape=thetaCurr[3], rate=thetaCurr[3], log=TRUE)))))
  }
  if(is.na(pacc[3])){
    pacc[3]=0
  }
  if(runif(1)<pacc[3]){
    thetaCurr[3]=thetaProp[3]
  } else {
    thetaProp[3]=thetaCurr[3]
  }
  
  #beta_e
  thetaProp[4]=runif(1, thetaCurr[4]-propLim[4], thetaCurr[4] + propLim[4])
  if(thetaProp[4]>0){
    pacc[4]=min(1, exp(logPriorDensityBetaE(thetaProp[4], betaEA, betaEB) +
                         sum(sapply(e[(tau+1):length(y)], dgamma, shape=thetaProp[4], rate=thetaProp[4], log=TRUE)) - (
                           logPriorDensityBetaE(thetaCurr[4], betaEA, betaEB) +
                             sum(sapply(e[(tau+1):length(y)], dgamma, shape=thetaCurr[4], rate=thetaCurr[4], log=TRUE)))))
  }
  
  if(is.na(pacc[4])){
    pacc[4]=0
  }
  if(runif(1)<pacc[4]){
    thetaCurr[4]=thetaProp[4]
  } else {
    thetaProp[4]=thetaCurr[4]
  }
  
  
  thetaProp[5]=rgamma(1, phiA + sum(y), rate=phiB + sum(RStates+SStates))
  thetaCurr[5]=thetaProp[5]
  pacc[5]=1
  
  return(list(thetaCurr, pacc))
}

######### set tuning parameters and initialization
N=50
rESSMin=0.5
its=50000
N0=N0
tau=5
state_lag=tau+1

RInit<-c(rep(0, 5), trunc(0.5*y[6:length(y)]))
SInit<-trunc(0.5*y)
SInit[1]=min(SInit[1], N0)
for(i in 2:length(y)){
  SInit[i]=min(SInit[i], SInit[i-1] + RInit[i-1])
}

thetaInit=numeric(5)
xx=SInit/c(N0, RInit[-length(y)] + SInit[-length(y)])
remove_xx=c(which(xx==0), which(is.na(xx)))
thetaInit[1]=mean(-log(xx[-remove_xx])/epsilon[-remove_xx])
thetaInit[2]=sum(RInit[-(1:tau)])/sum(c(N0, RInit[1:(length(y)-tau-1)] + SInit[1:(length(y)-tau-1)])*e[-(1:tau)]*exp(-c(N0, RInit[1:(length(y)-tau-1)] + SInit[1:(length(y)-tau-1)])/N0))
thetaInit[3]=1/var(epsilon)
thetaInit[4]=1/var(e[-(1:tau)])
thetaInit[5]=sum(y)/sum(RInit+SInit)


######### run algorithm 
theta<-matrix(rep(0, 5*its), nrow=its)
theta[1,]=thetaInit
RSamps<-matrix(rep(0, length(y)*its), nrow=its)
RSamps[1,]<-RInit
SSamps<-matrix(rep(0, length(y)*its), nrow=its)
SSamps[1,]<-SInit

timer<-system.time(
  for(i in 2:its){
    RSamps[i,]<-PGASR(N, rESSMin, y, theta[i-1,], e, epsilon, RSamps[i-1,], SSamps[i-1,])

    SSamps[i,]<-PGASS(N, rESSMin, y, theta[i-1,], e, epsilon, SSamps[i-1,], RSamps[i,])

    upTheta=thetaUpdate(SSamps[i,], RSamps[i,], theta[i-1,])
    theta[i,]=upTheta[[1]]
    
  }
)
