########### load packages 
library(invgamma)
library(truncnorm)

########### load function files 
source("../functions/run_MH.R")
source("../functions/gen_blocks.R")

########### simulate data 
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

########### define log observation and state density functions (fixed input and output)
logObsDensR<-function(y, states, t, theta){
  if(t<length(y)){
    p=exp(-theta[1]*epsilon[t+1])
    density=dpois(y[t], theta[5]*(S[t]+states[t]), log=TRUE) + dbinom(S[t+1], S[t]+states[t], p, log=TRUE)
  } else {
    density=dpois(y[t], theta[5]*(S[t]+states[t]), log=TRUE) 
  }
  density
  
} 
logObsDensS<-function(y, states, t, theta){
  if(t+tau+1<=length(y)){
    density=dpois(y[t], theta[5]*(R[t]+states[t]), log=TRUE) + dpois(R[t+tau+1], theta[2]*e[t+tau+1]*(R[t]+states[t])*exp(-(R[t]+states[t])/N0), log=T)
  } else {
    density=dpois(y[t], theta[5]*(R[t]+states[t]), log=TRUE)
    
  }
  density 
}

logStateDensR<-function(states, t, theta){
  if(t < stateLag){
    dens=log(1)
  } else {
    if(t==stateLag){
      dens=dpois(states[t], theta[2]*N0*e[t]*exp(-1), log=TRUE) 
    } else {
      if(t==length(y)){
        N=S[t-stateLag]+states[t-stateLag]
        dens=dpois(states[t], theta[2]*N*e[t]*exp(-N/N0), log=TRUE)
      } else {
        N=S[t-stateLag]+states[t-stateLag]
        dens=dpois(states[t], theta[2]*N*e[t]*exp(-N/N0), log=TRUE) 
      }
    }
  }
  dens
  
} 
logStateDensS<-function(states, t, theta){
  if(t==1){
    p=exp(-theta[1]*epsilon[t])
    dens=dbinom(states[t], N0, p, log=T) 
  } else {
    p=exp(-theta[1]*epsilon[t])
    dens=dbinom(states[t], states[t-1]+R[t-1], p, log=T)
  }
  
  dens
  
} 

########### define quantile and midpoint functions 
quantileFuncR<-function(perc, t, statesCurr, y, theta, varInfl, block){
  if(t<length(y)){
    lb=max(S[t+1]-S[t], 0)
    
  } else {
    lb=0
  }
  if(t<stateLag){
    quants=0
  } else {
    m=statesCurr[t]
    quants=pmax(qnorm(perc, m, sqrt(m*varInfl)), 0)
    #quants=unique(round(quants/2)*2)
    quants=round(quants)
    
  }
  quants=unique(pmax(lb, quants))
  sort(quants)
}
quantileFuncS<-function(perc, t, statesCurr, y, theta, varInfl, block){
  ub=numeric(length(na.omit(block)))
  lb=numeric(length(na.omit(block)))
  if(block[1]==1){
    ub[1]=N0
  } else {
    ub[1]=statesCurr[block[1]-1]+R[block[1]-1]
  }
  for(k in 2:length(ub)){
    ub[k]=ub[k-1]+R[block[k]-1]
  }
  if(block[length(na.omit(block))]==length(y)){
    lb[length(lb)]=0
  } else {
    lb[length(lb)]=max(0, statesCurr[block[length(na.omit(block))]+1]-R[block[length(na.omit(block))]])
  }
  
  for(k in (length(na.omit(block))-1):1){
    lb[k]=max(lb[k+1]-R[block[k]], 0)
  }
  
  m=statesCurr[t]
  pos=which(block==t)
  
  
  quants=pmax(qnorm(perc, m, sqrt(varInfl*m)), 0)

  quants=round(quants)
  quants=unique(pmin(quants, ub[pos]))
  quants=unique(pmax(lb[pos], quants))
  
  quants=sort(quants)
  quants 
  
}

midpointFuncR<-function(q, deltaE, t, statesCurr, block){
  if(t<length(y)){
    lb=max(S[t+1]-S[t], 0)
    
  } else {
    lb=0
  }
  if(t<stateLag){
    NAdapt=1
    mp=0
    lenCell=1
  } else {
    
    NAdapt<-length(q)+1
    
    mp<-numeric(NAdapt)
    lenCell<-numeric(NAdapt)
    
    mp[-c(1, NAdapt)]=(q[-1]+q[-(length(q))])/2
    lenCell[-c(1, NAdapt)]=c(q[-1]-q[-length(q)])
    
    if(NAdapt==2){
      lenCell[c(1, NAdapt)]=rep(2,2)
    } else {
      lenCell[c(1, NAdapt)]=pmax(round(mean(lenCell[-c(1, NAdapt)])/2)*2, 2)
    }
    mp[1]=pmax(as.integer(q[1]-lenCell[1]/2), 0)
    mp[NAdapt]<-as.integer(q[length(q)] + lenCell[NAdapt]/2)
    mp=ceiling(mp)
    
  }
  return(list("NAdapt"=NAdapt, "mp"=mp, "lenCell"=lenCell))
}
midpointFuncS<-function(q, deltaE, t, statesCurr, block){
  if(t<stateLag){
    NAdapt=1
    mp=0
    lenCell=1
  } else {
    NAdapt<-length(q)+1
    
    mp<-numeric(NAdapt)
    lenCell<-numeric(NAdapt)
    
    ub=numeric(length(na.omit(block)))
    if(block[1]==1){
      ub[1]=N0
    } else {
      ub[1]=statesCurr[block[1]-1]+R[block[1]-1]
    }
    for(k in 2:length(ub)){
      ub[k]=ub[k-1]+R[block[k]-1]
    }
    
    m=statesCurr[t]
    pos=which(block==t)

    mp[-c(1, NAdapt)]=(q[-1]+q[-(length(q))])/2
    lenCell[-c(1, NAdapt)]=c(q[-1]-q[-length(q)])
    
    if(NAdapt==2){
      lenCell[c(1, NAdapt)]=rep(2,2)
    } else {
      lenCell[c(1, NAdapt)]=pmax(round(mean(lenCell[-c(1, length(q))])/2)*2, 2)
    }
    mp[1]=pmax(as.integer(q[1]-lenCell[1]/2), 0)
    mp[NAdapt]<-min(as.integer(q[length(q)] + lenCell[NAdapt]/2), ub[pos])
    if(q[length(q)]==ub[pos]){
      mp[NAdapt]=ub[pos]+1
    }
    mp=ceiling(mp)
    
    
  }
  return(list("NAdapt"=NAdapt, "mp"=mp, "lenCell"=lenCell))
}

########### define midpoint integration functions 
midpointIntFuncSystemR<-function(mpoints, binLen, theta, t, stateLag){
  if(t<stateLag){
    dens=log(1)
  } else {
    if(t==stateLag){
      p=exp(-theta[1]*epsilon[t+1])
      dens=dbinom(S[t+1], S[t]+mpoints[[t]], p, log=TRUE) + dpois(mpoints[[t]], theta[2]*N0*e[t]*exp(-1), log=T) + log(binLen[[t]])
    } else {
      dens=matrix(rep(0, length(mpoints[[t-stateLag]])*length(mpoints[[t]])), ncol=length(mpoints[[t]]))
      for(k in 1:length(mpoints[[t-stateLag]])){
        N=S[t-stateLag]+mpoints[[t-stateLag]][k]
        if(t<length(y)){
          p=exp(-theta[1]*epsilon[t+1])
          dens[k,]=dpois(mpoints[[t]], theta[2]*N*e[t]*exp(-N/N0), log=TRUE) +
            dbinom(S[t+1], S[t]+mpoints[[t]], p, log=TRUE) + log(binLen[[t]]) + log(binLen[[t-stateLag]][k])
        } else {
          dens[k,]=dpois(mpoints[[t]], theta[2]*N*e[t]*exp(-N/N0), log=TRUE) + log(binLen[[t]]) + log(binLen[[t-stateLag]][k])
        }
      }
    }
  }
  
  
  dens
  
}
midpointIntFuncSystemS<-function(mpoints, binLen, theta, t, stateLag){
  if(t==1){
    dens=dbinom(mpoints[[t]], N0, exp(-theta[1]*epsilon[t]), log=T)
    if(length(which(mpoints[[t]]>N0))>0){
      dens[which(mpoints[[t]]>N0)]=-Inf 
    }
  } else {
    dens=matrix(rep(0, length(mpoints[[t-1]])*length(mpoints[[t]])), ncol=length(mpoints[[t]]))
    for(k in 1:length(mpoints[[t-1]])){
      N=R[t-1]+mpoints[[t-1]][k]
      p=exp(-theta[1]*epsilon[t])
      dens[k,]=dbinom(mpoints[[t]], N, p, log=T) + log(binLen[[t]]) + log(binLen[[t-1]][k])
      if(length(which(mpoints[[t]]>N))>0){
        dens[k,which(mpoints[[t]]>N)]=-Inf 
      }
    }
  }
  
  
  
  dens
  
}

midpointIntFuncObsR<-function(y, mpoints, binLen, theta, t){
  reformat<-function(mp, t) logObsDensR(y, c(numeric(t-1), mp), t, theta)
  sapply(mpoints[[t]], reformat, t=t) + log(binLen[[t]])
  
}
midpointIntFuncObsS<-function(y, mpoints, binLen, theta, t){
  reformat<-function(mp, t) logObsDensS(y, c(numeric(t-1), mp), t, theta)
  sapply(mpoints[[t]], reformat, t=t) + log(binLen[[t]])
  
}

correctThresh<-function(row, ind, thresh){
  if(length(which(row<thresh))>0){
    row[which(row<=thresh)]=thresh
  }
  row/sum(row)
}

########## proposal distributions for the states 
endGridCellPropR<-function(const, endPropVar, statesCurr, t, block){
  const+1+rpois(1, endPropVar)
}
endGridCellPropS<-function(const, endPropVar, statesCurr, t, block){
  const+1+rpois(1, endPropVar)
}
midGridCellProp<-function(lowerQuant, upperQuant){
  m=seq(lowerQuant+1, upperQuant, by=1)
  if(length(m)==1){
    m
  } else {
    sample(seq(lowerQuant+1, upperQuant, by=1), size=1)
  }
}
startGridCellPropR<-function(const, endPropVar, statesCurr, t, block){
  max(const-rpois(1, endPropVar), 0)
}
startGridCellPropS<-function(const, endPropVar, statesCurr, t, block){
  max(const-rpois(1, endPropVar), 0)
}

endGridCellDensR<-function(propState, const, endPropVar, statesCurr, t, block){
  dpois(propState-const-1, endPropVar, log=TRUE)
}
endGridCellDensS<-function(propState, const, endPropVar, statesCurr, t, block){
  dpois(propState-const-1, endPropVar, log=TRUE)
  }
midGridCellDens<-function(propState, lowerQuant, upperQuant){
  m=seq(lowerQuant+1, upperQuant, by=1)
  if(length(m)==1){
    log(1)
  } else {
    log(1/(upperQuant-lowerQuant))
  }
}
startGridCellDensR<-function(propState, const, endPropVar, statesCurr, t, block){
  if(propState==0){
    dens=log(1-ppois(const-1, endPropVar))
  } else {
    dens= dpois(const-propState, endPropVar, log=TRUE) 
  }
  dens
}
startGridCellDensS<-function(propState, const, endPropVar, statesCurr, t, block){
  if(propState==0){
    dens=log(1-ppois(const-1, endPropVar))
  } else {
    dens= dpois(const-propState, endPropVar, log=TRUE) 
  }
  dens
}

########## model parameter sampling
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
propEpsilon=0.2
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
N=50 #number of grid cells
lR=4 #block size in R
lS=4 #block size in S
q1R=0.1 #used to set span of quantile approaches (S)
qNR=1-q1R #used to set span of quantile approaches (S)
q1S=0.1 #used to set span of quantile approaches (S)
qNS=1-q1S #used to set span of quantile approaches (S)
varInflR=0.5 #used to set span of quantile approaches (R)
varInflS=0.5 #used to set span of quantile approaches (R)
thresh=0.01
endPropVarR=2 #proposal variance in extremal grid cells (R)
endPropVarS=2 #proposal variance in extremal grid cells (S)
approach=3 #specify approach to setting the grid cells 

its=50000

max.func<-function(el){
  max(0, el)
}
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
leny<-length(y)
bl<-genBlocks(lR, leny, tau+1)
blocksR=bl$blocks
bl<-genBlocks(lS, leny, 1)
blocksS=bl$blocks
tau=5

SStates<-matrix(rep(0, leny*its), nrow=its)
RStates<-matrix(rep(0, leny*its), nrow=its)
thetaCurr<-matrix(rep(0, length(thetaInit)*its), nrow=its)
SStates[1,]=SInit
RStates[1,]=RInit
thetaCurr[1,]<-thetaInit

timer<-system.time(
  for(i in 2:its){
    #update R states using PMPMH
    S=SStates[i-1,]
    logObsDens=logObsDensR
    logStateDens=logStateDensR
    quantileFunc<-quantileFuncR
    midpointFunc<-midpointFuncR
    endPropVar=endPropVarR
    midpointIntFuncSystem<-midpointIntFuncSystemR
    midpointIntFuncObs=midpointIntFuncObsR
    endGridCellProp=endGridCellPropR
    startGridCellProp=startGridCellPropR
    endGridCellDens=endGridCellDensR
    startGridCellDens=startGridCellDensR
    stateLag=tau+1
    varInfl=varInflR
    q1=q1R
    qN=qNR
    
    upStatesR<-PMPMH(RStates[i-1,], y, blocksR, N, q1, qN, varInfl, thetaCurr[i-1,], deltaE, thresh, approach, stateLag)
    RStates[i,]<-upStatesR$states
    
    #update S states 
    R=RStates[i,]
    logObsDens=logObsDensS
    logStateDens=logStateDensS
    quantileFunc<-quantileFuncS
    midpointFunc<-midpointFuncS
    endPropVar=endPropVarS
    midpointIntFuncSystem<-midpointIntFuncSystemS
    midpointIntFuncObs<-midpointIntFuncObsS
    endGridCellProp=endGridCellPropS
    startGridCellProp=startGridCellPropS
    endGridCellDens=endGridCellDensS
    startGridCellDens=startGridCellDensS
    stateLag=1
    varInfl=varInflS
    q1=q1S
    qN=qNS
   
    upStates<-PMPMH(SStates[i-1,], y, blocksS, N, q1, qN, varInfl, thetaCurr[i-1,], deltaE, thresh, approach, stateLag)
    SStates[i,]<-upStates$states
    
    #update parameters
    upTheta<-thetaUpdate(SStates[i,], RStates[i,], thetaCurr[i-1,])
    thetaCurr[i,]<-upTheta[[1]]

  }
)
