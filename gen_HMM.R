


norm_rows<-function(row){
  const<-max(row)
  if(const==-Inf){
    const==0
  }
  row<-exp(row-const)

  if(length(which(row==0))==length(row)){
    new.row=rep(1/length(row), length(row))
  } else {
    new.row=row/sum(row)
  }
 new.row
}

correct_thresh<-function(row, thresh){
  if(length(which(row<=thresh))>0){
  row[which(row<=thresh)]<-min(row[which(row>thresh)])
  }
  row/sum(row)
}

gen_HMM<-function(y, N, q1, qN, var.infl, states_curr, theta_curr, delta.e, thresh){
  len.y<-length(y)
  N.adapt<-numeric(len.y)
  perc<-seq(q1, qN, by=(qN-q1)/(N-2))

  q<-quantile_func(perc, 1, states_curr, y, theta_curr, var.infl)
  mps_extracted<-midpoint_func(q, delta.e, 1)

  N.adapt[1]<-mps_extracted$N.adapt
  mp<-mps_extracted$mp
  len.cell<-mps_extracted$len.cell

  quants=list(q)
  mpoints=list(mp)
  bin.len=list(len.cell)

  tr<-midpoint_int_func_system(mpoints, bin.len, theta_curr, 1)
  tr<-norm_rows(tr)
  tr<-correct_thresh(tr, thresh=thresh)

  ob<-midpoint_int_func_obs(y, mpoints, bin.len, theta_curr, 1)
  ob<-norm_rows(ob)
  ob<-correct_thresh(ob, thresh=thresh)

  transition=list(tr)
  observation=list(ob)

  for(t in 2:len.y){

    q=quantile_func(perc, t, states_curr, y, theta_curr, var.infl)
    mps_extracted<-midpoint_func(q, delta.e, t)

    N.adapt[t]<-mps_extracted$N.adapt
    mp<-mps_extracted$mp
    len.cell<-mps_extracted$len.cell

    quants[[t]]=q
    mpoints[[t]]=mp
    bin.len[[t]]=len.cell

    #use the other density function
    tr<-midpoint_int_func_system(mpoints, bin.len, theta_curr, t)
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


  return(list("Number bins"=N.adapt, "Transition"=transition, "Midpoints"=mpoints, "Quants"=quants, "Observation"=observation, "Bin length"=bin.len))
}

