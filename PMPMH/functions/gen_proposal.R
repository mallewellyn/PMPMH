findCurrIndex<-function(t, statesCurr, quants){
  findInterval(statesCurr[t], vec=c(-Inf, quants[[t]], Inf), left.open = TRUE)
}


sampleHMM<-function(statesCurr, y, block, grid, endPropVar, stateLag){
  transition=grid$transition
  observation=grid$observation
  NAdapt=grid$numberBins
  quants=grid$quants

  blockSize=length(na.omit(block))
  leny=length(y)

  indexCurr=numeric(leny)
  indexSamp<-numeric(leny)
  statesProp<-numeric(leny)
  filt<-as.list(NULL)

  indexCurr[na.omit(block)]=sapply(na.omit(block), findCurrIndex, statesCurr=statesCurr, quants=quants)


  filtTemp<-transition[[block[1]]]*observation[[block[1]]]
  filtTemp<-normRows(log(filtTemp))

  filt[[block[1]]]=filtTemp

  if(blockSize>1){
  for(i in 2:blockSize){
    filtTemp<-(filt[[block[i-1]]] %*% transition[[block[i]]])*observation[[block[i]]]
    filtTemp<-normRows(log(filtTemp))
    filt[[block[i]]]=filtTemp
  }
  }

  ## sampling
  if(block[blockSize]+stateLag>leny){
    smooth=filt[[block[blockSize]]]
  } else {
    smooth<-filt[[block[blockSize]]]*transition[[block[blockSize]+stateLag]]
  }
  smooth<-normRows(log(smooth))
  indexSamp[block[blockSize]]<-sample(1:NAdapt[block[blockSize]], size=1, prob=smooth)

  if(blockSize>1){
  for(i in rev(na.omit(block[-blockSize]))){
    smooth<-filt[[i]]*transition[[i+stateLag]][,indexSamp[i+stateLag]]
    smooth<-normRows(log(smooth))
    indexSamp[i]<-sample(1:NAdapt[i], size=1, prob=smooth)
  }
  }

  #proposal and current probability calculations
  probCalc<-function(t, index){
    log(transition[[t]][index[t-stateLag], index[t]]) + log(observation[[t]][index[t]])
  }

  logProbCurr<-log(transition[[block[1]]][indexCurr[block[1]]]) +
    log(observation[[block[1]]][indexCurr[block[1]]])
  if(blockSize>1){
    logProbCurr<-logProbCurr+sum(sapply(na.omit(block[-1]), probCalc, index=indexCurr))
  }
  if(block[blockSize]+stateLag<=leny){
    logProbCurr<-logProbCurr+log(transition[[block[blockSize]+stateLag]][indexCurr[block[blockSize]]])
  }

  logProbNew<-log(transition[[block[1]]][indexSamp[block[1]]]) +
    log(observation[[block[1]]][indexSamp[block[1]]])
  if(blockSize>1){
    logProbNew = logProbNew + sum(sapply(na.omit(block[-1]), probCalc, index=indexSamp))
  }

  if(block[blockSize]+stateLag<=leny){
    logProbNew = logProbNew+log(transition[[block[blockSize]+stateLag]][indexSamp[block[blockSize]]])
  }



  for(k in na.omit(block)){
    if(indexSamp[k]==NAdapt[k]){
      statesProp[k]<-endGridCellProp(quants[[k]][NAdapt[k]-1], endPropVar, statesCurr, k, block)
      logProbNew<-logProbNew + endGridCellDens(statesProp[k], quants[[k]][NAdapt[k]-1], endPropVar, statesCurr, k, block)
    }
    if(indexCurr[k]==NAdapt[k]){
      logProbCurr<-logProbCurr + endGridCellDens(statesCurr[k], quants[[k]][NAdapt[k]-1], endPropVar, statesCurr, k, block)
    }

    if(indexSamp[k]==1){
      statesProp[k]<-startGridCellProp(quants[[k]][1], endPropVar, statesCurr, k, block)
      logProbNew<-logProbNew + startGridCellDens(statesProp[k], quants[[k]][1], endPropVar, statesCurr, k, block)
    }
    if(indexCurr[k]==1){
      logProbCurr<-logProbCurr +startGridCellDens(statesCurr[k], quants[[k]][1], endPropVar, statesCurr, k, block)
    }

    if(indexSamp[k]!=1 & indexSamp[k]!=NAdapt[k]){
        statesProp[k]<-midGridCellProp(quants[[k]][indexSamp[k]-1], quants[[k]][indexSamp[k]])
      logProbNew<-logProbNew + midGridCellDens(statesProp[k], quants[[k]][indexSamp[k]-1], quants[[k]][indexSamp[k]])

    }
    if(indexCurr[k]!=1 & indexCurr[k]!=NAdapt[k]){
        logProbCurr<-logProbCurr + midGridCellDens(statesCurr[k], quants[[k]][indexCurr[k]-1], quants[[k]][indexCurr[k]])

    }

  }
  return(list("indexSamp"=indexSamp, "statesProp"=statesProp, "logProbNew"=logProbNew, "logProbCurr"=logProbCurr))
}


