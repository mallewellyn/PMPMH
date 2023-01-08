########### load packages 
library(invgamma)
library(truncnorm)

########### load function files 
source("../functions/run_MH.R")
source("../functions/gen_blocks.R")

########### simulate data 
set.seed(66)
n=600
x=numeric(n)
y=numeric(n)
priorMean=1
sigmaEpsilon=1
jumpProp=0.9
sigmaEtaSmall=1
sigmaEtaLarge=700

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


########### define log observation and state density functions (fixed input and output)
logObsDens<-function(y, states, t, theta){
  dnorm(y[t], states[t], sqrt(theta[4]), log=TRUE)
} 
logStateDens<-function(states, t, theta){
  if(t==1){
    log(theta[1]*dnorm(states[t], 1, sqrt(theta[2])) + (1-theta[1])*dnorm(states[t], 1, sqrt(theta[3])))
  } else {
    log(theta[1]*dnorm(states[t], states[t-1], sqrt(theta[2])) + (1-theta[1])*dnorm(states[t], states[t-1], sqrt(theta[3])))
  }
} 

########### define quantile and midpoint functions 
quantileFunc<-function(perc, t, statesCurr, y, theta, varInfl, block){
  quants=qnorm(perc, statesCurr[t], sqrt(varInfl))
  unique(quants)
}
midpointFunc<-function(q, deltaE, t, statesCurr, block){
  NAdapt<-length(q)+1
  
  mp<-numeric(NAdapt)
  lenCell<-numeric(NAdapt)
  
  mp[-c(1, NAdapt)]=(q[-1]+q[-(length(q))])/2
  lenCell[-c(1, NAdapt)]=q[-1]-q[-length(q)]
  
  avLenCell=mean(lenCell[-c(1, NAdapt)])
  
  mp[1]<-q[1]-avLenCell/2
  mp[NAdapt]<-q[length(q)] + avLenCell/2
  lenCell[c(1, NAdapt)]<-avLenCell
  
  return(list("NAdapt"=NAdapt, "mp"=mp, "lenCell"=lenCell))
}

########### define midpoint integration functions 
midpointIntFuncSystem<-function(mpoints, binLen, theta, t, stateLag){
  reformat<-function(mpMinus, mp) logStateDens(c(numeric(t-2), mpMinus, mp), t, theta)
  vReformat<-Vectorize(reformat, vectorize.args = c("mpMinus", "mp"))
  
  if(t==1){
    dens=sapply(mpoints[[1]], logStateDens, t=1, theta=theta) + log(binLen[[1]])
  } else {
    dens=outer(mpoints[[t-stateLag]], mpoints[[t]], vReformat)
    dens=dens+binLen[[t-stateLag]]
    dens=t(t(dens) + binLen[[t]])
  }
  
  
  dens
  
}
midpointIntFuncObs<-function(y, mpoints, binLen, theta, t){
  reformat<-function(mp, t) logObsDens(y, c(numeric(t-1), mp), t, theta)
  
  sapply(mpoints[[t]], reformat, t=t) + log(binLen[[t]])
  
}

correctThresh<-function(row, ind, thresh){
  if(length(which(row<thresh))>0){
    row[which(row<=thresh)]=thresh
  }
  row/sum(row)
}
########## proposal distributions for the states 
endGridCellProp<-function(const, endPropVar, statesCurr, t, block){
  rtruncnorm(1, a=const, b=Inf, mean=const, sqrt(endPropVar))
}
startGridCellProp<-function(const, endPropVar, statesCurr, t, block){
  rtruncnorm(1, a=-Inf, b=const, mean=const, sqrt(endPropVar))
}
midGridCellProp<-function(lowerQuant, upperQuant){
  runif(1, lowerQuant, upperQuant)
}

endGridCellDens<-function(propState, const, endPropVar, statesCurr, t, block){
  log(dtruncnorm(propState, a=const, b=Inf, mean=const, sqrt(endPropVar)))
}
startGridCellDens<-function(propState, const, endPropVar, statesCurr, t, block){
  log(dtruncnorm(propState, a=-Inf, b=const, mean=const, sqrt(endPropVar)))
}
midGridCellDens<-function(propState, lowerQuant, upperQuant){
  dunif(propState, lowerQuant, upperQuant, log=TRUE)
}

########## model parameter sampling
betaSigma=2
gammaSigma=2
betaSmall=2
gammaSmall=2
betaLarge=2
gammaLarge=700
propLim=c(0.3, 2, 160, NA) #p, small, large, epsilon

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
  
  return(list("upTheta"=thetaCurr, "propTheta"=thetaProp, "Pacc"=pacc))
  
}


######### set tuning parameters and initialization
N=10 #number of grid cells
l=4 #block size
q1=0.1 #used to set span of quantile approaches
qN=1-q1 #used to set span of quantile approaches
varInfl=1.4 #used to set span of quantile approaches
thresh=0.01 #threshold for correction
endPropVar=5 #proposal variance in extremal grid cells 
approach=3 #specify approach to setting the grid cells 


its=10000
stateLag=1
xInit=y
thetaInit=c(0.5, 5, 5, 5)

######### run algorithm
leny<-length(y)
bl<-genBlocks(l, leny, stateLag)
blocks=bl$blocks

statesCurr<-matrix(rep(0, leny*its), nrow=its)
thetaCurr<-matrix(rep(0, length(thetaInit)*its), nrow=its)
statesCurr[1,]<-xInit
thetaCurr[1,]<-thetaInit

timer<-system.time(
  for(i in 2:its){
    #update x samples using PMPMH
    upStates<-PMPMH(statesCurr[i-1,], y, blocks, N, q1, qN, varInfl, thetaCurr[i-1,], deltaE, thresh, approach, stateLag)
    statesCurr[i,]<-upStates$states
    
    upTheta<-thetaUpdate(statesCurr[i,], thetaCurr[i-1,], propLim)
    thetaCurr[i,]<-upTheta$upTheta
  }
)

