findCurrIndex<-function(t, statesCurr, quants){
  findInterval(statesCurr[t], vec=c(-Inf, quants[[t]], Inf), left.open = TRUE)
}

normRows<-function(row){
  const=max(row)
  if(const==Inf){
    row[which(row==Inf)]=0
    row=exp(row)
  }
  if(const==-Inf){
    row=rep(1, length(row))
  }
  if(const != Inf & const !=-Inf){
    row=exp(row-const)
  }
 row/sum(row)
}


genHMM<-function(block, y, N, q1, qN, varInfl, statesCurr, thetaCurr, deltaE, thresh, stateLag){
  quants=as.list(NULL)
  mpoints=as.list(NULL)
  binLen=as.list(NULL)
  transition=as.list(NULL)
  observation=as.list(NULL)

  if(block[1]-stateLag>0){
    mpoints[[block[1]-stateLag]]=statesCurr[block[1]-stateLag]
    binLen[[block[1]-stateLag]]=1
  }

  lenBlock<-length(block)
  leny<-length(y)
  NAdapt<-numeric(leny)
  perc<-seq(q1, qN, by=(qN-q1)/(N-2))

  q<-quantileFunc(perc, block[1], statesCurr, y, thetaCurr, varInfl, na.omit(block))
  mpsExtracted<-midpointFunc(q, deltaE, block[1], statesCurr, na.omit(block))

  NAdapt[block[1]]<-mpsExtracted$NAdapt
  mp<-mpsExtracted$mp
  lenCell<-mpsExtracted$lenCell

  quants[[block[1]]]=q
  mpoints[[block[1]]]=mp
  binLen[[block[1]]]=lenCell

  ###
  ind=findCurrIndex(block[1], statesCurr, quants)

  tr=midpointIntFuncSystem(mpoints, binLen, thetaCurr, block[1], stateLag)
  tr<-normRows(tr)
  tr<-correctThresh(tr, ind, thresh)

  ob<-midpointIntFuncObs(y, mpoints, binLen, thetaCurr, block[1])
  ob<-normRows(ob)
  ob<-correctThresh(ob, ind, thresh)

  transition[[block[1]]]=tr
  observation[[block[1]]]=ob

  t=block[1]

  if(length(na.omit(block[-1]))>0){
  for(t in na.omit(block[-1])){

    if(t<max(na.omit(block))){
    q=quantileFunc(perc, t, statesCurr, y, thetaCurr, varInfl, na.omit(block))
    mpsExtracted<-midpointFunc(q, deltaE, t, statesCurr, na.omit(block))
    } else {
      q=quantileFunc(perc, t, statesCurr, y, thetaCurr, varInfl, na.omit(block))
      mpsExtracted<-midpointFunc(q, deltaE, t, statesCurr, na.omit(block))
    }

    NAdapt[t]<-mpsExtracted$NAdapt
    mp<-mpsExtracted$mp
    lenCell<-mpsExtracted$lenCell

    quants[[t]]=q
    mpoints[[t]]=mp
    binLen[[t]]=lenCell

    #use the other density function
    ind=findCurrIndex(t, statesCurr, quants)

    tr<-midpointIntFuncSystem(mpoints, binLen, thetaCurr, t, stateLag)
    if(is.matrix(tr)){
      tr<-t(apply(tr, MARGIN=1, normRows))
      tr=t(apply(tr, MARGIN=1, correctThresh, ind=ind, thresh=thresh))
    } else {
      tr<-normRows(tr)
      tr<-correctThresh(tr, ind, thresh)
    }


    ob<-midpointIntFuncObs(y, mpoints, binLen, thetaCurr, t)
    ob<-normRows(ob)
    ob<-correctThresh(ob, ind, thresh)

    transition[[t]]=tr
    observation[[t]]=ob

  }
  }

  if(t+stateLag<=length(y)){
    i=block[length(block)] + stateLag

    mpoints[[i]]=statesCurr[i]
    binLen[[i]]=1

    tr<-midpointIntFuncSystem(mpoints, binLen, thetaCurr, i, stateLag)
    if(is.matrix(tr)){
         tr<-t(apply(tr, MARGIN=1, normRows))
    } else {
      tr<-normRows(tr)
    }

    transition[[i]]=tr

  }


  return(list("numberBins"=NAdapt, "transition"=transition, "midpoints"=mpoints, "quants"=quants, "observation"=observation, "binLength"=binLen))
}

