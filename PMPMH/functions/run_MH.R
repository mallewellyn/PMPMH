source("../functions/gen_HMM.R")
source("../functions/gen_proposal.R")

genPacc<-function(block, prop, statesCurr, thetaCurr, stateLag){
  statesProp<-statesCurr
  statesProp[block]<-prop$statesProp[block]
  
  if(block[1]<2*stateLag){
    tEvalStates<-c(block, max(block)+stateLag)
  }
  if(max(block)+stateLag>length(y)){
    tEvalStates<-block
  }
  if(min(block)!=stateLag && max(block)+stateLag<=length(y)){
    tEvalStates<-c(block, max(block)+stateLag)
  }
  
  pacc<-min(1, exp(prop$logProbCurr + sum(sapply(tEvalStates, logStateDens, states=statesProp, theta=thetaCurr)) +
                     sum(sapply(block, logObsDens, y=y, states=statesProp, theta=thetaCurr)) -
                     (prop$logProbNew + sum(sapply(tEvalStates, logStateDens, states=statesCurr, theta=thetaCurr)) +
                        sum(sapply(block, logObsDens, y=y, states=statesCurr, theta=thetaCurr)))))
  
  if(is.na(pacc)){
    pacc=0
  }
  pacc
}
PMPMH<-function(statesCurr, y, blocks, N, q1, qN, varInfl, thetaCurr, deltaE, thresh, approach, stateLag){
  
  leny=length(y)
  
  paccStates<-numeric(dim(blocks)[1])
  propsCurr=paccStates
  propsNew=paccStates
  statesProp<-numeric(leny)
  NCells<-numeric(dim(blocks)[1])
  
  if(approach!=3){
    grid=genHMM(seq(1, leny, by=1), y, N, q1, qN, varInfl, statesCurr, thetaCurr, deltaE, thresh, stateLag)
  }
  
  for(l in 1:dim(blocks)[1]){
    if(approach==3){
      grid=genHMM(na.omit(blocks[l,]), y, N, q1, qN, varInfl, statesCurr, thetaCurr, deltaE, thresh, stateLag)
      NCells[l]=mean(grid$numberBins[na.omit(blocks[l,])])
    }
    proposal<-sampleHMM(statesCurr, y, na.omit(blocks[l,]), grid, endPropVar, stateLag)
    paccStates[l]<-genPacc(na.omit(blocks[l,]), proposal, statesCurr, thetaCurr, stateLag)
    statesProp[na.omit(blocks[l,])]=proposal$statesProp[na.omit(blocks[l,])]
    if(runif(1)<paccStates[l]){
      statesCurr[na.omit(blocks[l,])]<-proposal$statesProp[na.omit(blocks[l,])]
    }
  }
  
  return(list("states"=statesCurr, "proposedStates"=statesProp, "statesPacc"=paccStates, "NCells"=NCells, "propProbCurr"=propsCurr, "propProbNew"=propsNew))
}