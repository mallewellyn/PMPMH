


norm_rows<-function(row){
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

correct_thresh<-function(row, thresh){
  if(length(which(row<thresh))>0){
  row[which(row<=thresh)]<-min(row[which(row>thresh)])
  }
  row/sum(row)
}

gen_HMM<-function(block, y, N, q1, qN, var.infl, states_curr, theta_curr, delta.e, thresh, state_lag){
  quants=as.list(NULL)
  mpoints=as.list(NULL)
  bin.len=as.list(NULL)
  transition=as.list(NULL)
  observation=as.list(NULL)

  len.block<-length(block)
  len.y<-length(y)
  N.adapt<-numeric(len.y)
  perc<-seq(q1, qN, by=(qN-q1)/(N-2))

  q<-quantile_func(perc, block[1], states_curr, y, theta_curr, var.infl)
  mps_extracted<-midpoint_func(q, delta.e, 1)

  N.adapt[block[1]]<-mps_extracted$N.adapt
  mp<-mps_extracted$mp
  len.cell<-mps_extracted$len.cell

  quants[[block[1]]]=q
  mpoints[[block[1]]]=mp
  bin.len[[block[1]]]=len.cell

  if(block[1]>1){
    mpoints[[block[1]-state_lag]]=states_curr[block[1]-state_lag]
    bin.len[[block[1]-state_lag]]=1
  }

  tr=midpoint_int_func_system(mpoints, bin.len, theta_curr, block[1], state_lag)
  tr<-norm_rows(tr)
  tr<-correct_thresh(tr, thresh=thresh)

  ob<-midpoint_int_func_obs(y, mpoints, bin.len, theta_curr, block[1])
  ob<-norm_rows(ob)
  ob<-correct_thresh(ob, thresh=thresh)

  transition[[block[1]]]=tr
  observation[[block[1]]]=ob

  for(t in block[-1]){

    q=quantile_func(perc, t, states_curr, y, theta_curr, var.infl)
    mps_extracted<-midpoint_func(q, delta.e, t)

    N.adapt[t]<-mps_extracted$N.adapt
    mp<-mps_extracted$mp
    len.cell<-mps_extracted$len.cell

    quants[[t]]=q
    mpoints[[t]]=mp
    bin.len[[t]]=len.cell

    #use the other density function
    tr<-midpoint_int_func_system(mpoints, bin.len, theta_curr, t, state_lag)
    if(is.matrix(tr)){
      tr<-t(apply(tr, MARGIN=1, norm_rows))
      tr<-t(apply(tr, MARGIN=1, correct_thresh, thresh=thresh))
    } else {
      tr<-norm_rows(tr)
      tr<-correct_thresh(tr, thresh=thresh)
    }

    ob<-midpoint_int_func_obs(y, mpoints, bin.len, theta_curr, t)
    ob<-norm_rows(ob)
    ob<-correct_thresh(ob, thresh=thresh)

    transition[[t]]=tr
    observation[[t]]=ob

  }

  if(t<length(y)){
    mpoints[[block[length(block)] + state_lag]]=states_curr[block[length(block)] + state_lag]
    bin.len[[block[length(block)] + state_lag]]=1

    tr<-midpoint_int_func_system(mpoints, bin.len, theta_curr, block[length(block)]+state_lag, state_lag)
    if(is.matrix(tr)){
      tr<-t(apply(tr, MARGIN=1, norm_rows))
      tr<-t(apply(tr, MARGIN=1, correct_thresh, thresh=thresh))
    } else {
      tr<-norm_rows(tr)
      tr<-correct_thresh(tr, thresh=thresh)
    }

    transition[[block[length(block)] + state_lag]]=tr
  }


  return(list("Number bins"=N.adapt, "Transition"=transition, "Midpoints"=mpoints, "Quants"=quants, "Observation"=observation, "Bin length"=bin.len))
}

