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
  
  if(approach==1){
    grid_temp=genHMM(seq(1, length(y), by=1), y, N, q1, qN, varInfl, statesCurr, thetaCurr, deltaE, thresh, stateLag)
  }
  
  for(l in 1:dim(blocks)[1]){
    
    if(approach==1){
      grid=grid_temp
      if(blocks[l,1]-stateLag>0){
        grid$midpoints[[blocks[l,1]-stateLag]]=statesCurr[blocks[l,1]-stateLag]
        grid$binLength[[blocks[l,1]-stateLag]]=1
        
        ind=findCurrIndex(blocks[l,1]-stateLag, statesCurr, grid$quants)
        grid$transition[[blocks[l,1]]]=grid$transition[[blocks[l,1]]][ind,]
      }
      
      bl=na.omit(blocks[l,])
      if(bl[length(bl)]+stateLag<=length(y)){
        grid$midpoints[[bl[length(bl)]+stateLag]]=statesCurr[bl[length(bl)]+stateLag]
        grid$binLength[[bl[length(bl)]+stateLag]]=1
        grid$transition[[bl[length(bl)]+stateLag]]=t(as.matrix(rep(1, grid$numberBins[[bl[length(bl)]]]), ncol=1))
        }
      
    } else {
      grid=genHMM(na.omit(blocks[l,]), y, N, q1, qN, varInfl, statesCurr, thetaCurr, deltaE, thresh, stateLag)
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