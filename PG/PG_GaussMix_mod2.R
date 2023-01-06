########### load packages 
library(invgamma)

######### simulate data 
set.seed(100)
n=1000
x=numeric(n)
y=numeric(n)
priorMean=1
sigmaEpsilon=10
jumpProp=0.99
sigmaEtaSmall=1
sigmaEtaLarge=10000

u=runif(n)
if(u[1]<jumpProp){
  x[1]=rnorm(1, priorMean, sqrt(sigmaEtaSmall))
} else {
  x[1]=rnorm(1, priorMean, sqrt(sigmaEtaLarge))
}
y[1]=rnorm(1, x[1], sqrt(sigmaEpsilon))
for(i in 2:n){
  if(u[i]<jumpProp){
    x[i]=rnorm(1, x[i-1], sqrt(sigmaEtaSmall))
  } else {
    x[i]=rnorm(1, x[i-1], sqrt(sigmaEtaLarge))
  }
  y[i]=rnorm(1, x[i], sqrt(sigmaEpsilon))
}
set.seed(NULL)

######### define random draws from prior (state distribution), and evaluations of observation, prior densities
priorRandInit<-function(t, theta){
  if(runif(1)<theta[1]){
    rnorm(1, priorMean, sqrt(theta[2]))
  } else {
    rnorm(1, priorMean, sqrt(theta[3]))
  }
}
priorRandOtherwise<-function(xtMinus1, theta){
  if(runif(1)<theta[1]){
    rnorm(1, xtMinus1, sqrt(theta[2]))
  } else {
    rnorm(1, xtMinus1, sqrt(theta[3]))
  }
}

logLikelihood<-function(t, xt, theta){
  dnorm(y[t], xt, sqrt(theta[4]), log=TRUE)
}
logPrior<-function(child, theta, partMinus){
  log(theta[1]*dnorm(child, partMinus, sqrt(theta[2])) + (1-theta[1])*dnorm(child, partMinus, sqrt(theta[3])))
}

########### PG function
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

PG<-function(N, rESSMin, y, theta, xInit){
  
  T=length(y)
  particles<-matrix(rep(0, N*T), ncol=T)
  weights<-matrix(rep(0, N*T), ncol=T)
  A<-matrix(rep(0, N*T), ncol=T)
  normWeights<-matrix(rep(0, N*T), ncol=T)
  
  #time 1
  particles[2:N,1]<-sapply(rep(1, N-1), priorRandInit, theta=theta)
  particles[1,1]=xInit[1]
  
  logWeights<-logLikelihood(1, particles[,1], theta)
  if(length(which(logWeights==-Inf))==N){
    logWeights=rep(log(1), N)
  }
  const=max(logWeights)
  weights[,1]=exp(logWeights-const)
  normWeights[,1]=weights[,1]/sum(weights[,1])
  
  for(t in 2:T){
    A[1,t]=1
    if(((1/sum(normWeights[,t-1]^2))/N)<rESSMin){
      A[2:N,t]<-multinomialResample(normWeights[,t-1])
      weights[,t-1]=1
    } else {
      A[2:N,t]=2:N
    }
    
    particles[2:N,t]=sapply(particles[,t-1][A[-1,t]], priorRandOtherwise, theta=theta)
    particles[1,t]=xInit[t]
    
    particles[,1:t]=cbind(particles[A[,t], 1:(t-1)], particles[,t])
    ##include weights[A[,t],t-1]?
    logWeights<-log(weights[A[,t],t-1])+logLikelihood(t, particles[,t], theta)
    if(length(which(logWeights==-Inf))==N){
      logWeights=rep(log(1), N)
    }
    const=max(logWeights)
    weights[,t]<-exp(logWeights-const)
    normWeights[,t]<-weights[,t]/sum(weights[,t])
    
    
  }
  
  xx<-sample(multinomialResample(normWeights[,T]), size=1)
  samp<-particles[xx, ]
  
  
  samp
}

########## model parameter sampling
betaSigma=2
gammaSigma=10
betaSmall=2
gammaSmall=2
betaLarge=2
gammaLarge=10000
propLim=c(0.02, 0.5, 20000, NA) #p, small, large, epsilon

thetaUpdate<-function(states, theta, propLim){
  thetaCurr=theta
  thetaProp=theta
  pacc=numeric(length(thetaCurr))
  
  #p
  thetaProp[1]=runif(1, thetaCurr[1]-propLim[1]/2, thetaCurr[1]+propLim[1]/2)
  
  if(thetaProp[1]>0 && thetaProp[1]<1){
    pacc[1]=min(1, exp(log(thetaProp[1]*dnorm(states[1], priorMean, sqrt(thetaCurr[2])) + (1-thetaProp[1])*dnorm(states[1], priorMean, sqrt(thetaCurr[3]))) + 
                         sum(log(thetaProp[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[2])) + (1-thetaProp[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[3]))))
                       - (log(thetaCurr[1]*dnorm(states[1], priorMean, sqrt(thetaCurr[2])) + (1-thetaCurr[1])*dnorm(states[1], priorMean, sqrt(thetaCurr[3]))) + 
                            sum(log(thetaCurr[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[2])) + (1-thetaCurr[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[3])))))))
    
    if(is.na(pacc[1])){
      pacc[1]=0
    }
    if(runif(1)<pacc[1]){
      thetaCurr[1]=thetaProp[1]
    } else {
      thetaProp[1]=thetaCurr[1]
    }
  }
  
  #sigma_eta small
  thetaProp[2]=runif(1, thetaCurr[2]-propLim[2]/2, thetaCurr[2]+propLim[2]/2)
  
  if(thetaProp[2]>0 && thetaProp[2]<thetaCurr[3]){
    
    pacc[2]=min(1, exp(dinvgamma(thetaProp[2], shape=betaSmall, rate=gammaSmall, log=T) + 
                         log(thetaCurr[1]*dnorm(states[1], priorMean, sqrt(thetaProp[2])) + (1-thetaCurr[1])*dnorm(states[1], priorMean, sqrt(thetaCurr[3]))) +
                         sum(log(thetaCurr[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaProp[2])) +(1-thetaCurr[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[3])))) -
                         (dinvgamma(thetaCurr[2], shape=betaSmall, rate=gammaSmall, log=T) + 
                            log(thetaCurr[1]*dnorm(states[1], priorMean, sqrt(thetaCurr[2])) + (1-thetaCurr[1])*dnorm(states[1], priorMean, sqrt(thetaCurr[3]))) +
                            sum(log(thetaCurr[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[2])) +(1-thetaCurr[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[3]))))))) 
    
    
    if(is.na(pacc[2])){
      pacc[2]=0
    }
    if(runif(1)<pacc[2]){
      thetaCurr[2]=thetaProp[2]
    } else {
      thetaProp[2]=thetaCurr[2]
    }
  }
  
  #sigma_eta large 
  thetaProp[3]=runif(1, thetaCurr[3]-propLim[3]/2, thetaCurr[3]+propLim[3]/2)
  
  if(thetaProp[3]>0){
    pacc[3]=min(1, exp(dinvgamma(thetaProp[3], shape=betaLarge, rate=gammaLarge, log=T) + 
                         log(thetaCurr[1]*dnorm(states[1], priorMean, sqrt(thetaCurr[2])) + (1-thetaCurr[1])*dnorm(states[1], priorMean, sqrt(thetaProp[3]))) +
                         sum(log(thetaCurr[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[2])) +(1-thetaCurr[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaProp[3])))) -
                         (dinvgamma(thetaCurr[3], shape=betaLarge, rate=gammaLarge, log=T) + 
                            log(thetaCurr[1]*dnorm(states[1], priorMean, sqrt(thetaCurr[2])) + (1-thetaCurr[1])*dnorm(states[1], priorMean, sqrt(thetaCurr[3]))) +
                            sum(log(thetaCurr[1]*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[2])) +(1-thetaCurr[1])*dnorm(states[-1], states[-length(states)], sqrt(thetaCurr[3]))))))) 
    
    if(is.na(pacc[3])){
      pacc[3]=0
    }
    if(runif(1)<pacc[3]){
      thetaCurr[3]=thetaProp[3]
    } else {
      thetaProp[3]=thetaCurr[3]
    }
    
  }
  
  #sigma_epsilon
  thetaCurr[4]=rinvgamma(1, shape=betaSigma+length(y)/2, rate=gammaSigma + (1/2)*sum((y-states)^2))
  thetaProp[4]=thetaCurr[4]
  
  return(list("upTheta"=thetaCurr, "propTheta"=thetaProp, "pacc"=pacc))
  
}

######### set tuning parameters and initialization
N=10
rESSMin=0.5
its=10000
xInit=y
thetaInit=c(0.5, 5, 5, 5)



########### run algorithm
xSamps<-matrix(rep(0, length(y)*its), nrow=its)
theta<-matrix(rep(0, 4*its), nrow=its)
xSamps[1,]=xInit
theta[1,]=thetaInit

timer<-system.time(
    for(i in 2:its){
      #update x samples using PG
      xSamps[i,]<-PG(N, rESSMin, y, theta[i-1,], xSamps[i-1,])
      
      #update parameters using M-H
      upTheta=thetaUpdate(xSamps[i,], theta[i-1,], propLim)
      theta[i,]=upTheta$upTheta
}
)
