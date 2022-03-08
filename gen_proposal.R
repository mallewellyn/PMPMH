
sample_HMM<-function(states_curr, y, block, grid, end.prop.var, state_lag){
  transition=grid$Transition
  observation=grid$Observation
  N.adapt=grid$`Number bins`
  quants=grid$Quants

  blocksize=length(block)
  len.y=length(y)

  index.curr<-numeric(len.y)
  for(k in block){
    index.curr[k]<-find_curr_index(states_curr[k], quants[[k]])
  }

  index.samp<-numeric(len.y)
  states_prop<-numeric(len.y)
  filt<-as.list(NULL)

  filt_temp<-transition[[block[1]]]*observation[[block[1]]]
  filt_temp<-norm_rows(filt_temp)

  filt[[block[1]]]=filt_temp

  for(i in 2:blocksize){
    filt_temp<-(filt[[block[i-1]]] %*% transition[[block[i]]])*observation[[block[i]]]
    filt_temp<-norm_rows(filt_temp)
    filt[[block[i]]]=filt_temp
  }

  ## sampling
  if(len.y==block[blocksize]){
    smooth=filt[[block[blocksize]]]
  } else {
    smooth<-filt[[block[blocksize]]]*transition[[block[blocksize]+state_lag]]
  }
  smooth<-norm_rows(smooth)
  index.samp[block[blocksize]]<-which((runif(1)-(cumsum(smooth)/cumsum(smooth)[length(smooth)])<0))[1]

  for(i in rev(block[-blocksize])){
    smooth<-filt[[i]]*transition[[i+state_lag]][,index.samp[i+state_lag]]
    smooth<-norm_rows(smooth)
    index.samp[i]<-which((runif(1)-(cumsum(smooth)/cumsum(smooth)[length(smooth)])<0))[1]
  }

  #proposal and current probability calculations
  prob_calc<-function(t, index){
    log(transition[[t]][index[t-state_lag], index[t]]) + log(observation[[t]][index[t]])
  }

  log.prob.curr<-log(transition[[block[1]]][index.curr[block[1]]]) +
    log(observation[[block[1]]][index.curr[block[1]]]) +
    sum(sapply(block[-1], prob_calc, index=index.curr))
  if(len.y!=block[blocksize]){
    log.prob.curr<-log.prob.curr+log(transition[[block[blocksize]+state_lag]][index.curr[block[blocksize]]])
  }

  log.prob.new<-log(transition[[block[1]]][index.samp[block[1]]]) +
    log(observation[[block[1]]][index.samp[block[1]]]) +
    sum(sapply(block[-1], prob_calc, index=index.samp))

  if(len.y!=block[blocksize]){
  log.prob.new<-log.prob.new+log(transition[[block[blocksize]+state_lag]][index.samp[block[blocksize]]])
  }



  for(k in block){
    if(index.samp[k]==N.adapt[k]){
      states_prop[k]<-end_grid_cell_prop(quants[[k]][N.adapt[k]-1], end.prop.var)
      log.prob.new<-log.prob.new + end_grid_cell_dens(states_prop[k], quants[[k]][N.adapt[k]-1], end.prop.var)
    }
    if(index.curr[k]==N.adapt[k]){
      log.prob.curr<-log.prob.curr + end_grid_cell_dens(states_curr[k], quants[[k]][N.adapt[k]-1], end.prop.var)
    }

    if(index.samp[k]==1){
      states_prop[k]<-start_grid_cell_prop(quants[[k]][1], end.prop.var)
      log.prob.new<-log.prob.new + start_grid_cell_dens(states_prop[k], quants[[k]][1], end.prop.var)
    }
    if(index.curr[k]==1){
      log.prob.curr<-log.prob.curr +start_grid_cell_dens(states_curr[k], quants[[k]][1], end.prop.var)
    }

    if(index.samp[k]!=1 & index.samp[k]!=N.adapt[k]){
        states_prop[k]<-mid_grid_cell_prop(quants[[k]][index.samp[k]-1], quants[[k]][index.samp[k]])
      log.prob.new<-log.prob.new+mid_grid_cell_dens(states_prop[k], quants[[k]][index.samp[k]-1], quants[[k]][index.samp[k]])

    }
    if(index.curr[k]!=1 & index.curr[k]!=N.adapt[k]){
        log.prob.curr<-log.prob.curr + mid_grid_cell_dens(states_curr[k], quants[[k]][index.curr[k]-1], quants[[k]][index.curr[k]])

    }

  }
  return(list("index.samp"=index.samp, "states_prop"=states_prop, "log.prob.new"=log.prob.new, "log.prob.curr"=log.prob.curr))
}

